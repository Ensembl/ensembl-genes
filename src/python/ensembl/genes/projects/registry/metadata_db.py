"""
Data fetcher for the Ensembl Metadata database.
"""
import pymysql
import logging
from typing import Optional, Dict, Any
from ensembl.genes.projects.models import GenomeMetadata

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class MetadataDbClient:
    """Queries the Ensembl Metadata Schema to build GenomeMetadata."""
    
    def __init__(self, host: str, port: int, user: str, dbname: str):
        self.host = host
        self.port = port
        self.user = user
        self.dbname = dbname

    def fetch_by_identifier(self, identifier: str) -> Optional[GenomeMetadata]:
        """
        Identifier can be a Genome UUID, a core DB name, or an Assembly Accession.
        """
        # Determine the field to search
        if "-" in identifier and len(identifier) == 36:
            where_clause = "genome.genome_uuid = %s"
        elif "core" in identifier or "_cm_" in identifier:
            where_clause = "dataset_source.name = %s"
        else:
            where_clause = "assembly.accession = %s"

        query = f"""
            SELECT 
                genome.genome_uuid,
                dataset_source.name AS dbname,
                assembly.accession,
                organism.scientific_name,
                organism.strain,
                organism.taxonomy_id,
                assembly.name AS assembly_name
            FROM genome
            JOIN assembly ON genome.assembly_id = assembly.assembly_id
            JOIN organism ON genome.organism_id = organism.organism_id
            JOIN genome_dataset ON genome.genome_id = genome_dataset.genome_id
            JOIN dataset ON genome_dataset.dataset_id = dataset.dataset_id
            JOIN dataset_source ON dataset.dataset_source_id = dataset_source.dataset_source_id
            WHERE dataset.name = 'genebuild'
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
            logger.error(f"MySQL error querying metadata DB: {e}")
            return None
        finally:
            if 'conn' in locals() and conn.open:
                conn.close()

        if not row:
            logger.warning(f"No metadata found for identifier: {identifier}")
            return None
            
        species_name = row['scientific_name'].replace(" ", "_").lower() if row['scientific_name'] else ""
        
        # Proper alternate pairing: look for assemblies of the exact same taxonomy & strain
        # but with a different accession. Returns the rapid-release URL.
        alternate_url = None
        if row.get('strain'):
            alt_query = """
                SELECT organism.url_name
                FROM genome
                JOIN organism ON genome.organism_id = organism.organism_id
                JOIN assembly ON genome.assembly_id = assembly.assembly_id
                WHERE organism.taxonomy_id = %s
                  AND organism.strain = %s
                  AND assembly.accession != %s
                LIMIT 1
            """
            try:
                conn = pymysql.connect(host=self.host, user=self.user, port=self.port, database=self.dbname)
                with conn.cursor(pymysql.cursors.DictCursor) as cursor:
                    cursor.execute(alt_query, (row['taxonomy_id'], row['strain'], row['accession']))
                    alt_row = cursor.fetchone()
                    if alt_row and alt_row.get('url_name'):
                        alternate_url = f"https://rapid.ensembl.org/{alt_row['url_name']}/Info/Index"
            except Exception as e:
                logger.error(f"Failed to lookup alternate haplotype: {e}")
            finally:
                if 'conn' in locals() and conn.open:
                    conn.close()

        return GenomeMetadata(
            genome_uuid=row['genome_uuid'],
            dbname=row['dbname'],
            accession=row['accession'],
            species_name=species_name,
            assembly_name=row['assembly_name'],
            strain=row['strain'],
            taxon_id=row['taxonomy_id'],
            alternate_of=alternate_url
        )
