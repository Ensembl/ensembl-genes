"""
Data fetcher for the Genebuild Pre-release Tracker database (gb_schema).
Replaces fragile 'information_schema' matching for tracking genomes in active production.
"""
import pymysql
import logging
from typing import Optional, List
from ensembl.genes.projects.models import GenomeMetadata
from ensembl.genes.projects.config import ProjectConfig

logger = logging.getLogger(__name__)

class GbTrackerClient:
    """Queries the new gb_schema DB to build GenomeMetadata for pre-release cores."""
    
    def __init__(self, host: str, port: int, user: str, dbname: str = "gb_assembly_metadata"):
        self.host = host
        self.port = port
        self.user = user
        self.dbname = dbname

    def fetch_by_identifier(self, identifier: str) -> Optional[GenomeMetadata]:
        """
        Looks up a pre-release core database name or accession in the tracking tables.
        """
        if "_core_" in identifier or "_cm_" in identifier:
            separator = "_core_" if "_core_" in identifier else "_cm_"
            species_search = identifier.split(separator)[0]
            where_clause = "REPLACE(LOWER(s.scientific_name), ' ', '_') = %s"
            query_param = species_search
        else:
            where_clause = "gs.gca_accession = %s"
            query_param = identifier

        query = f"""
            SELECT 
                NULL AS genome_uuid,
                NULL AS core_dbname,
                gs.gca_accession AS assembly_accession,
                s.scientific_name AS species_name,
                s.species_taxon_id AS taxon_id,
                a.asm_name AS assembly_name,
                gs.annotation_method
            FROM genebuild_status gs
            JOIN assembly a ON gs.assembly_id = a.assembly_id
            JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
            WHERE gs.gb_status IN ('in_progress', 'pre_released')
              AND {where_clause}
            LIMIT 1
        """
        
        try:
            conn = pymysql.connect(
                host=self.host, user=self.user, port=self.port, database=self.dbname
            )
            with conn.cursor(pymysql.cursors.DictCursor) as cursor:
                cursor.execute(query, (query_param,))
                row = cursor.fetchone()
        except pymysql.Error as e:
            logger.error(f"MySQL error querying gb_schema: {e}")
            return None
        finally:
            if 'conn' in locals() and conn.open:
                conn.close()

        if not row:
            return None

        # Pre-release DBs don't have beta links yet, but they have FTPs constructed locally
        raw_method = row.get('annotation_method')
        method_map = {
            'braker': 'BRAKER2',
            'full_genebuild': 'Ensembl Genebuild',
            'mixed_strategy_build': 'Mixed Strategy Build',
        }
        anno_method = method_map.get(raw_method) if raw_method else raw_method

        return GenomeMetadata(
            genome_uuid=row.get('genome_uuid') or "unknown",
            dbname=row['core_dbname'] or "unknown",
            accession=row['assembly_accession'],
            species_name=row['species_name'],
            assembly_name=row['assembly_name'],
            taxon_id=row.get('taxon_id'),
            annotation_method=anno_method,
            is_on_rapid=False,
            is_on_beta=False,
            is_released=False
        )

    def fetch_project_pre_releases(self, config: ProjectConfig) -> List[GenomeMetadata]:
        """
        Discovers pre-release genomes from the GB registry scoping exclusively by the explicit
        ProjectConfig rules (e.g. bioproject_scoping or custom_group_scoping).
        """
        if not config.bioproject_scoping and not config.custom_group_scoping:
            logger.info(f"Project '{config.project_name}' has no defined pre-release scoping bounds in GB registry.")
            return []

        where_clauses = []
        params = []
        joins = ""

        if config.bioproject_scoping:
            joins += " JOIN bioproject bp ON gs.assembly_id = bp.assembly_id"
            joins += " JOIN main_bioproject mb ON bp.bioproject_id = mb.bioproject_id"
            
            in_clause = ', '.join(['%s'] * len(config.bioproject_scoping))
            where_clauses.append(f"mb.bioproject_name IN ({in_clause})")
            params.extend(config.bioproject_scoping)
            
        elif config.custom_group_scoping:
            joins += " JOIN custom_group cg ON gs.assembly_id = cg.item AND cg.group_type = 'assembly'"
            
            in_clause = ', '.join(['%s'] * len(config.custom_group_scoping))
            where_clauses.append(f"cg.group_name IN ({in_clause})")
            params.extend(config.custom_group_scoping)

        where_cond = " OR ".join(where_clauses) if len(where_clauses) > 1 else where_clauses[0]

        query = f"""
            SELECT 
                NULL AS genome_uuid,
                NULL AS core_dbname,
                gs.gca_accession AS assembly_accession,
                s.scientific_name AS species_name,
                s.species_taxon_id AS taxon_id,
                a.asm_name AS assembly_name,
                gs.annotation_method,
                MAX(CASE WHEN am.metrics_name = 'genebuild.busco' THEN am.metrics_value END) AS busco_score,
                MAX(CASE WHEN am.metrics_name = 'genebuild.busco_dataset' THEN am.metrics_value END) AS busco_lineage
            FROM genebuild_status gs
            JOIN assembly a ON gs.assembly_id = a.assembly_id
            JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
            {joins}
            LEFT JOIN annotation_metrics am ON gs.assembly_id = am.assembly_id AND gs.genebuild_status_id = am.genebuild_status_id
            WHERE gs.gb_status = 'pre_released'
              AND ({where_cond})
            GROUP BY gs.gca_accession, s.scientific_name, s.species_taxon_id, a.asm_name, gs.annotation_method
        """

        try:
            conn = pymysql.connect(
                host=self.host, user=self.user, port=self.port, database=self.dbname
            )
            with conn.cursor(pymysql.cursors.DictCursor) as cursor:
                cursor.execute(query, tuple(params))
                rows = cursor.fetchall()
        except pymysql.Error as e:
            logger.error(f"MySQL error querying gb_schema for project pre-releases: {e}")
            return []
        finally:
            if 'conn' in locals() and conn.open:
                conn.close()

        pre_releases = []
        method_map = {
            'braker': 'BRAKER2',
            'full_genebuild': 'Ensembl Genebuild',
            'mixed_strategy_build': 'Mixed Strategy Build',
        }

        for row in rows:
            raw_method = row.get('annotation_method')
            anno_method = method_map.get(raw_method) if raw_method else raw_method

            pre_releases.append(
                GenomeMetadata(
                    genome_uuid=row.get('genome_uuid') or "unknown",
                    dbname=row.get('core_dbname') or "unknown",
                    accession=row['assembly_accession'],
                    species_name=row['species_name'],
                    assembly_name=row['assembly_name'],
                    taxon_id=row.get('taxon_id'),
                    annotation_method=anno_method,
                    busco_score=row.get('busco_score'),
                    busco_lineage=row.get('busco_lineage'),
                    is_on_rapid=False,
                    is_on_beta=False,
                    is_released=False
                )
            )

        return pre_releases
