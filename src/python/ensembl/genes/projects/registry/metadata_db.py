"""
Data fetcher for the Ensembl Metadata database.
"""

import logging
from typing import Dict, List, Optional

import pymysql

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

    def _fetch_taxonomy_lineage(self, taxonomy_id: int) -> List[str]:
        """Fetch taxonomy classification names for a given taxonomy_id.

        Queries the ``organism_classification`` table in the metadata DB,
        which stores the NCBI lineage as (organism_id, cls_group_id, rank)
        rows linked back to the organism table.

        Falls back to an empty list if the table does not exist or the
        query fails — this is expected on some DB versions.

        Returns names ordered leaf → root (most specific first).
        """
        if not taxonomy_id:
            return []

        # Strategy: try several known schema shapes in order of likelihood.
        # The metadata DB schema varies across Ensembl releases; we try
        # the most common query patterns and silently fall back on failure.
        queries = [
            # New metadata schema: taxonomy_node parent traversal
            # Walk from the organism's taxonomy_id up through parent nodes.
            """
            WITH RECURSIVE lineage AS (
                SELECT taxonomy_id, name, parent_id
                FROM taxonomy_node
                WHERE taxonomy_id = %s
                UNION ALL
                SELECT tn.taxonomy_id, tn.name, tn.parent_id
                FROM taxonomy_node tn
                JOIN lineage l ON tn.taxonomy_id = l.parent_id
                WHERE tn.taxonomy_id != tn.parent_id
            )
            SELECT name FROM lineage
            WHERE taxonomy_id != %s
            """,
            # Alternate: flat organism_classification table
            """
            SELECT oc.name
            FROM organism o
            JOIN organism_classification oc ON o.organism_id = oc.organism_id
            WHERE o.taxonomy_id = %s
            ORDER BY oc.classification_order ASC
            """,
        ]

        for query in queries:
            try:
                conn = pymysql.connect(
                    host=self.host,
                    user=self.user,
                    port=self.port,
                    database=self.dbname,
                )
                with conn.cursor() as cursor:
                    if query.count("%s") == 2:
                        cursor.execute(query, (taxonomy_id, taxonomy_id))
                    else:
                        cursor.execute(query, (taxonomy_id,))
                    rows = cursor.fetchall()
                conn.close()

                if rows:
                    names = [r[0] for r in rows if r[0]]
                    logger.debug(
                        "Fetched %d lineage entries for taxon_id=%d from metadata DB",
                        len(names),
                        taxonomy_id,
                    )
                    return names
            except pymysql.Error as e:
                # Table doesn't exist or query is incompatible — try next
                logger.debug(
                    "Taxonomy lineage query failed for taxon_id=%d: %s",
                    taxonomy_id,
                    e,
                )
                try:
                    if "conn" in locals() and conn.open:
                        conn.close()
                except Exception:  # pylint: disable=broad-exception-caught
                    pass
                continue

        return []

    def get_genome_uuids_by_accessions(self, accessions: List[str]) -> Dict[str, str]:
        """Look up genome UUIDs for a list of assembly accessions.

        Parameters
        ----------
        accessions:
            GCA/GCF accession strings to look up.

        Returns
        -------
        Dict mapping accession → genome_uuid for accessions found in the
        metadata DB.  Missing accessions are silently omitted.
        """
        if not accessions:
            return {}

        # Deduplicate and filter empties
        unique = sorted(set(a for a in accessions if a))
        if not unique:
            return {}

        placeholders = ", ".join(["%s"] * len(unique))
        query = f"""
            SELECT DISTINCT
                assembly.accession,
                genome.genome_uuid
            FROM genome
            JOIN assembly ON genome.assembly_id = assembly.assembly_id
            JOIN genome_dataset ON genome.genome_id = genome_dataset.genome_id
            JOIN dataset ON genome_dataset.dataset_id = dataset.dataset_id
            WHERE dataset.name = 'genebuild'
              AND assembly.accession IN ({placeholders})
        """

        result: Dict[str, str] = {}
        try:
            conn = pymysql.connect(
                host=self.host,
                user=self.user,
                port=self.port,
                database=self.dbname,
            )
            with conn.cursor() as cursor:
                cursor.execute(query, tuple(unique))
                for row in cursor.fetchall():
                    acc, guuid = row[0], row[1]
                    if acc and guuid:
                        result[acc] = guuid
        except pymysql.Error as e:
            logger.error("MySQL error looking up genome UUIDs: %s", e)
        finally:
            if "conn" in locals() and conn.open:
                conn.close()

        logger.debug(
            "Resolved %d / %d accessions to genome UUIDs from metadata DB.",
            len(result),
            len(unique),
        )
        return result

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
                assembly.name AS assembly_name,
                MAX(CASE WHEN attribute.name = 'genebuild.method_display' THEN dataset_attribute.value END) AS annotation_method,
                MAX(CASE WHEN attribute.name = 'genebuild.annotation_source' THEN dataset_attribute.value END) AS annotation_source,
                MAX(CASE WHEN attribute.name = 'genebuild.last_geneset_update' THEN dataset_attribute.value END) AS geneset_date,
                MAX(CASE WHEN attribute.name = 'genebuild.busco' THEN dataset_attribute.value END) AS busco_score,
                MAX(CASE WHEN attribute.name = 'genebuild.busco_dataset' THEN dataset_attribute.value END) AS busco_lineage
            FROM genome
            JOIN assembly ON genome.assembly_id = assembly.assembly_id
            JOIN organism ON genome.organism_id = organism.organism_id
            JOIN genome_dataset ON genome.genome_id = genome_dataset.genome_id
            JOIN dataset ON genome_dataset.dataset_id = dataset.dataset_id
            JOIN dataset_source ON dataset.dataset_source_id = dataset_source.dataset_source_id
            LEFT JOIN dataset_attribute ON dataset.dataset_id = dataset_attribute.dataset_id
            LEFT JOIN attribute ON dataset_attribute.attribute_id = attribute.attribute_id 
                AND attribute.name IN ('genebuild.method_display', 'genebuild.annotation_source', 'genebuild.last_geneset_update', 'genebuild.busco', 'genebuild.busco_dataset')
            WHERE dataset.name = 'genebuild'
              AND {where_clause}
            GROUP BY genome.genome_uuid, dataset_source.name, assembly.accession, organism.scientific_name, organism.strain, organism.taxonomy_id, assembly.name
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
            logger.error("MySQL error querying metadata DB: %s", e)
            return None
        finally:
            if "conn" in locals() and conn.open:
                conn.close()

        if not row:
            logger.warning("No metadata found for identifier: %s", identifier)
            return None

        species_name = row["scientific_name"] if row["scientific_name"] else ""

        # Proper alternate pairing: look for assemblies of the exact same taxonomy & strain
        # but with a different accession. Returns the rapid-release URL.
        alternate_url = None
        if row.get("strain"):
            alt_query = """
                SELECT genome.url_name
                FROM genome
                JOIN organism ON genome.organism_id = organism.organism_id
                JOIN assembly ON genome.assembly_id = assembly.assembly_id
                WHERE organism.taxonomy_id = %s
                  AND organism.strain = %s
                  AND assembly.accession != %s
                LIMIT 1
            """
            try:
                conn = pymysql.connect(
                    host=self.host, user=self.user, port=self.port, database=self.dbname
                )
                with conn.cursor(pymysql.cursors.DictCursor) as cursor:
                    cursor.execute(
                        alt_query, (row["taxonomy_id"], row["strain"], row["accession"])
                    )
                    alt_row = cursor.fetchone()
                    if alt_row and alt_row.get("url_name"):
                        alternate_url = (
                            "https://rapid.ensembl.org/"
                            f"{alt_row['url_name']}/Info/Index"
                        )
            except Exception as e:  # pylint: disable=broad-exception-caught
                logger.error("Failed to lookup alternate haplotype: %s", e)
            finally:
                if "conn" in locals() and conn.open:
                    conn.close()

        # Attempt to fetch taxonomy lineage from the metadata DB
        taxonomy_lineage = self._fetch_taxonomy_lineage(row["taxonomy_id"])

        return GenomeMetadata(
            genome_uuid=row["genome_uuid"],
            dbname=row["dbname"],
            accession=row["accession"],
            species_name=species_name,
            assembly_name=row["assembly_name"],
            strain=row["strain"],
            taxon_id=row["taxonomy_id"],
            taxonomy_lineage=taxonomy_lineage or None,
            alternate_of=alternate_url,
            annotation_source=row.get("annotation_source"),
            annotation_method=row.get("annotation_method"),
            annotation_date=row.get("geneset_date"),
            busco_score=row.get("busco_score"),
            busco_lineage=row.get("busco_lineage"),
            is_released=True,
        )
