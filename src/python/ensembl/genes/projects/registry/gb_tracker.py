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
    
    def __init__(self, host: str, port: int, user: str, dbname: str = "gb_schema"):
        self.host = host
        self.port = port
        self.user = user
        self.dbname = dbname

    def fetch_by_identifier(self, identifier: str) -> Optional[GenomeMetadata]:
        """
        Looks up a pre-release core database name or accession in the tracking tables.
        """
        if "core" in identifier or "_cm_" in identifier:
            # `gb_schema` does not natively store the final core_dbname.
            # We match by accession primarily. If a core DB name is passed, we might fail unless we parse it.
            # For safety, we search the `gca_accession` or `scientific_name`.
            where_clause = "gs.gca_accession = %s"
        else:
            where_clause = "gs.gca_accession = %s"

        query = f"""
            SELECT 
                NULL AS genome_uuid,
                NULL AS core_dbname,
                gs.gca_accession AS assembly_accession,
                s.scientific_name AS species_name,
                a.asm_name AS assembly_name
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
                cursor.execute(query, (identifier,))
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
        return GenomeMetadata(
            genome_uuid=row.get('genome_uuid') or "unknown",
            dbname=row['core_dbname'],
            accession=row['assembly_accession'],
            species_name=row['species_name'],
            assembly_name=row['assembly_name'],
            is_on_rapid=False,
            is_on_beta=False
        )
