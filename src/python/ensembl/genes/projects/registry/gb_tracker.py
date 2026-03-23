"""
Data fetcher for the Genebuild Pre-release Tracker database (gb_schema).
Replaces fragile 'information_schema' matching for tracking genomes in active production.
"""
import pymysql
import logging
from typing import Optional
from ensembl.genes.projects.models import GenomeMetadata

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
            annotation_method=anno_method,
            is_on_rapid=False,
            is_on_beta=False,
            is_released=False
        )
