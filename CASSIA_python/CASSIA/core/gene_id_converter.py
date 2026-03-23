"""
Gene ID Converter for CASSIA.

This module provides automatic conversion of Ensembl and Entrez gene IDs
to gene symbols. Uses the mygene library for ID mapping, with Ensembl REST
API as a fallback.

Features:
- Auto-detect species from Ensembl ID prefix (ENSG=human, ENSMUSG=mouse, etc.)
- Strip version suffixes from Ensembl IDs (e.g., ENSG00000141510.17)
- Warn when ID prefix doesn't match the specified species
- Fall back to Ensembl REST API when mygene fails
"""

import re
import warnings
from typing import List, Tuple, Dict, Optional, Set

# Gene ID patterns - also captures optional version suffix
ENSEMBL_PATTERN = re.compile(r'^ENS[A-Z]*G\d{11}(\.\d+)?$', re.IGNORECASE)
ENSEMBL_STRIP_VERSION = re.compile(r'^(ENS[A-Z]*G\d{11})\.\d+$', re.IGNORECASE)
ENTREZ_PATTERN = re.compile(r'^\d+$')

# Ensembl ID prefix -> (species name for mygene, common name)
ENSEMBL_PREFIX_TO_SPECIES = {
    'ENSG':      ('human', 'human'),
    'ENSMUSG':   ('mouse', 'mouse'),
    'ENSRNOG':   ('rat', 'rat'),
    'ENSSSCG':   ('pig', 'pig'),
    'ENSBTAG':   ('bovine', 'cow'),
    'ENSCAFG':   ('dog', 'dog'),
    'ENSFELG':   ('cat', 'cat'),           # Felis catus (not always in mygene)
    'ENSGALG':   ('chicken', 'chicken'),    # Gallus gallus
    'ENSDARG':   ('zebrafish', 'zebrafish'),
    'ENSXETG':   ('frog', 'frog'),          # Xenopus tropicalis
    'ENSOARG':   ('sheep', 'sheep'),        # Ovis aries
    'ENSECAG':   ('horse', 'horse'),        # Equus caballus
    'ENSOCUG':   ('rabbit', 'rabbit'),      # Oryctolagus cuniculus
    'ENSMMUG':   ('macaque', 'macaque'),     # Macaca mulatta (rhesus)
    'ENSPTRG':   ('chimpanzee', 'chimpanzee'),
    'ENSGGOG':   ('gorilla', 'gorilla'),
    'ENSPANG':   ('olive baboon', 'baboon'),
    'FBgn':      ('fruitfly', 'fruit fly'),  # Drosophila
    'WBGene':    ('nematode', 'C. elegans'),
}

# Try to import mygene
MYGENE_AVAILABLE = False
try:
    import mygene
    MYGENE_AVAILABLE = True
except ImportError:
    mygene = None


def strip_ensembl_version(gene_id: str) -> str:
    """Strip version suffix from Ensembl ID (e.g., ENSG00000141510.17 -> ENSG00000141510)."""
    m = ENSEMBL_STRIP_VERSION.match(gene_id)
    return m.group(1) if m else gene_id


def detect_species_from_ids(ensembl_ids: List[str]) -> Optional[str]:
    """
    Auto-detect species from Ensembl ID prefixes.

    Looks at the prefix pattern (e.g., ENSG, ENSMUSG, ENSSSCG) and maps
    to a mygene-compatible species name.

    Returns:
        Species string if detected, None if ambiguous or unknown.
    """
    if not ensembl_ids:
        return None

    # Sort prefixes longest-first to avoid e.g. ENSG matching before ENSGALG
    sorted_prefixes = sorted(ENSEMBL_PREFIX_TO_SPECIES.items(), key=lambda x: -len(x[0]))

    detected = {}
    for eid in ensembl_ids:
        eid_upper = eid.upper()
        for prefix, (species, _) in sorted_prefixes:
            if eid_upper.startswith(prefix.upper()):
                detected[species] = detected.get(species, 0) + 1
                break

    if not detected:
        return None

    # Return the most common species
    best = max(detected, key=detected.get)
    return best


def warn_species_mismatch(user_species: str, detected_species: str, ensembl_ids: List[str]):
    """Warn user if specified species doesn't match Ensembl ID prefix."""
    if user_species.lower() == detected_species.lower():
        return

    # Get the common name for display
    common_name = detected_species
    for _, (sp, cn) in ENSEMBL_PREFIX_TO_SPECIES.items():
        if sp == detected_species:
            common_name = cn
            break

    sample = ensembl_ids[:3]
    warnings.warn(
        f"[CASSIA] Species mismatch: you specified species='{user_species}', "
        f"but Ensembl ID prefixes (e.g., {', '.join(sample)}) suggest '{common_name}'. "
        f"Auto-correcting to species='{detected_species}' for better conversion results. "
        f"If this is wrong, pass the correct species parameter explicitly.",
        UserWarning
    )


def _ensembl_rest_convert(ensembl_ids: List[str], species: str = 'human') -> Dict[str, str]:
    """
    Fallback: convert Ensembl IDs to gene symbols via Ensembl REST API.

    Uses POST /lookup/id to batch-query Ensembl IDs.
    Returns a mapping of ensembl_id -> symbol for successful lookups.
    """
    try:
        import urllib.request
        import json
    except ImportError:
        return {}

    if not ensembl_ids:
        return {}

    mapping = {}
    # Ensembl REST API accepts up to 1000 IDs per request
    batch_size = 1000
    for i in range(0, len(ensembl_ids), batch_size):
        batch = ensembl_ids[i:i + batch_size]
        try:
            url = "https://rest.ensembl.org/lookup/id"
            payload = json.dumps({"ids": batch}).encode('utf-8')
            req = urllib.request.Request(
                url,
                data=payload,
                headers={
                    "Content-Type": "application/json",
                    "Accept": "application/json"
                },
                method="POST"
            )
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read().decode('utf-8'))

            for eid, info in data.items():
                if isinstance(info, dict) and info.get('display_name'):
                    mapping[eid] = info['display_name']
        except Exception:
            # Silently continue — this is a fallback
            continue

    return mapping


def detect_id_type(gene_id: str) -> str:
    """
    Detect the type of gene identifier.

    Args:
        gene_id: A gene identifier string

    Returns:
        str: One of 'ensembl', 'entrez', or 'symbol'
    """
    if ENSEMBL_PATTERN.match(gene_id):
        return 'ensembl'
    elif ENTREZ_PATTERN.match(gene_id) and len(gene_id) > 2:
        return 'entrez'
    else:
        return 'symbol'


def classify_markers(markers: List[str]) -> Dict[str, List[str]]:
    """
    Classify a list of markers by their ID type.

    Args:
        markers: List of gene identifiers

    Returns:
        dict: Dictionary with keys 'ensembl', 'entrez', 'symbol'
              and lists of corresponding IDs
    """
    classified = {
        'ensembl': [],
        'entrez': [],
        'symbol': []
    }

    for marker in markers:
        id_type = detect_id_type(marker)
        classified[id_type].append(marker)

    return classified


def convert_ensembl_to_symbols(
    ensembl_ids: List[str],
    species: str = 'human'
) -> Tuple[Dict[str, str], List[str]]:
    """
    Convert Ensembl gene IDs to gene symbols using mygene, with Ensembl REST fallback.

    Handles version-suffixed IDs (e.g., ENSG00000141510.17) by stripping versions
    before querying and mapping results back to original IDs.

    Args:
        ensembl_ids: List of Ensembl gene IDs (may include version suffixes)
        species: Species for gene lookup ('human', 'mouse', 'pig', etc.)

    Returns:
        Tuple of (mapping dict {original_id: symbol}, list of failed original IDs)
    """
    if not ensembl_ids:
        return {}, []

    # Step 1: Strip version suffixes and build original -> stripped mapping
    original_to_stripped = {}
    stripped_to_originals = {}  # stripped can map to multiple originals
    for eid in ensembl_ids:
        stripped = strip_ensembl_version(eid)
        original_to_stripped[eid] = stripped
        if stripped not in stripped_to_originals:
            stripped_to_originals[stripped] = []
        stripped_to_originals[stripped].append(eid)

    unique_stripped = list(stripped_to_originals.keys())

    # Step 2: Auto-detect species from ID prefixes
    detected_species = detect_species_from_ids(unique_stripped)
    effective_species = species
    if detected_species and detected_species.lower() != species.lower():
        warn_species_mismatch(species, detected_species, unique_stripped)
        effective_species = detected_species

    # Step 3: Try mygene first
    stripped_mapping = {}
    mygene_failed = list(unique_stripped)  # assume all failed initially

    if MYGENE_AVAILABLE:
        try:
            mg = mygene.MyGeneInfo()
            results = mg.querymany(
                unique_stripped,
                scopes='ensembl.gene',
                fields='symbol',
                species=effective_species,
                returnall=True
            )

            mygene_failed = []
            for result in results.get('out', []):
                query_id = result.get('query', '')
                symbol = result.get('symbol')
                if symbol:
                    stripped_mapping[query_id] = symbol
                else:
                    mygene_failed.append(query_id)

            # Add any IDs that weren't in the results
            queried_ids = {r.get('query', '') for r in results.get('out', [])}
            for sid in unique_stripped:
                if sid not in queried_ids:
                    mygene_failed.append(sid)

        except Exception as e:
            warnings.warn(
                f"mygene query failed: {str(e)}. Trying Ensembl REST API fallback.",
                UserWarning
            )
            mygene_failed = unique_stripped

    # Step 4: Ensembl REST API fallback for IDs that mygene couldn't convert
    if mygene_failed:
        rest_mapping = _ensembl_rest_convert(mygene_failed, effective_species)
        stripped_mapping.update(rest_mapping)
        still_failed = [sid for sid in mygene_failed if sid not in rest_mapping]
        if rest_mapping:
            n_rest = len(rest_mapping)
            n_remaining = len(still_failed)
            print(f"[CASSIA] Ensembl REST API recovered {n_rest} additional gene symbol(s)"
                  f"{f' ({n_remaining} still unresolved)' if n_remaining > 0 else ''}.")
    else:
        still_failed = []

    # Step 5: Map back to original IDs (with version suffixes)
    mapping = {}
    failed = []
    for eid in ensembl_ids:
        stripped = original_to_stripped[eid]
        if stripped in stripped_mapping:
            mapping[eid] = stripped_mapping[stripped]
        else:
            failed.append(eid)

    return mapping, failed


def convert_entrez_to_symbols(
    entrez_ids: List[str],
    species: str = 'human'
) -> Tuple[Dict[str, str], List[str]]:
    """
    Convert Entrez/NCBI gene IDs to gene symbols using mygene.

    Args:
        entrez_ids: List of Entrez gene IDs (as strings)
        species: Species for gene lookup ('human' or 'mouse')

    Returns:
        Tuple of (mapping dict, list of failed IDs)
    """
    if not MYGENE_AVAILABLE:
        return {}, entrez_ids

    if not entrez_ids:
        return {}, []

    try:
        mg = mygene.MyGeneInfo()

        # Query mygene for the Entrez IDs
        results = mg.querymany(
            entrez_ids,
            scopes='entrezgene',
            fields='symbol',
            species=species,
            returnall=True
        )

        mapping = {}
        failed = []

        for result in results.get('out', []):
            query_id = result.get('query', '')
            symbol = result.get('symbol')

            if symbol:
                mapping[query_id] = symbol
            else:
                failed.append(query_id)

        # Add any IDs that weren't in the results
        queried_ids = {r.get('query', '') for r in results.get('out', [])}
        for eid in entrez_ids:
            if eid not in queried_ids:
                failed.append(eid)

        return mapping, failed

    except Exception as e:
        warnings.warn(
            f"Error converting Entrez IDs: {str(e)}. "
            "Returning original IDs.",
            UserWarning
        )
        return {}, entrez_ids


def convert_gene_ids(
    markers: List[str],
    species: str = 'human',
    auto_convert: bool = True
) -> Tuple[List[str], Dict[str, any]]:
    """
    Automatically detect and convert Ensembl/Entrez IDs to gene symbols.

    This is the main entry point for gene ID conversion.

    Args:
        markers: List of gene identifiers (may be mixed types)
        species: Species for gene lookup ('human' or 'mouse')
        auto_convert: Whether to attempt automatic conversion

    Returns:
        Tuple of:
            - List of converted markers (gene symbols where possible)
            - Conversion info dict with keys:
                - 'converted_count': Number of IDs successfully converted
                - 'failed_count': Number of IDs that couldn't be converted
                - 'ensembl_converted': List of (original, symbol) tuples
                - 'entrez_converted': List of (original, symbol) tuples
                - 'failed_ids': List of IDs that couldn't be converted
                - 'conversion_performed': Whether any conversion was done
    """
    info = {
        'converted_count': 0,
        'failed_count': 0,
        'ensembl_converted': [],
        'entrez_converted': [],
        'failed_ids': [],
        'conversion_performed': False
    }

    if not auto_convert:
        return markers, info

    # Classify markers by type
    classified = classify_markers(markers)

    # If no Ensembl or Entrez IDs, return as-is
    if not classified['ensembl'] and not classified['entrez']:
        return markers, info

    # Warn if mygene not available (Ensembl REST fallback still works for Ensembl IDs)
    if not MYGENE_AVAILABLE:
        if classified['entrez']:
            warnings.warn(
                f"Detected {len(classified['entrez'])} Entrez IDs but 'mygene' package is not installed. "
                "Entrez conversion requires mygene. Install with: pip install mygene\n"
                "Ensembl IDs will still be converted via Ensembl REST API.",
                UserWarning
            )
        if not classified['ensembl']:
            # Only Entrez IDs and no mygene — can't do anything
            return markers, info

    info['conversion_performed'] = True

    # Convert Ensembl IDs (handles version stripping, species detection, REST fallback internally)
    ensembl_mapping, ensembl_failed = convert_ensembl_to_symbols(
        classified['ensembl'], species
    )

    # Convert Entrez IDs
    entrez_mapping, entrez_failed = convert_entrez_to_symbols(
        classified['entrez'], species
    )

    # Build converted marker list preserving order
    converted_markers = []
    for marker in markers:
        id_type = detect_id_type(marker)

        if id_type == 'ensembl' and marker in ensembl_mapping:
            symbol = ensembl_mapping[marker]
            converted_markers.append(symbol)
            info['ensembl_converted'].append((marker, symbol))
            info['converted_count'] += 1
        elif id_type == 'entrez' and marker in entrez_mapping:
            symbol = entrez_mapping[marker]
            converted_markers.append(symbol)
            info['entrez_converted'].append((marker, symbol))
            info['converted_count'] += 1
        elif id_type in ('ensembl', 'entrez'):
            # Failed to convert - keep original
            converted_markers.append(marker)
            info['failed_ids'].append(marker)
            info['failed_count'] += 1
        else:
            # Already a symbol
            converted_markers.append(marker)

    # Log conversion summary
    if info['converted_count'] > 0:
        msg_parts = []
        if info['ensembl_converted']:
            msg_parts.append(f"{len(info['ensembl_converted'])} Ensembl IDs")
        if info['entrez_converted']:
            msg_parts.append(f"{len(info['entrez_converted'])} Entrez IDs")

        conversion_msg = (
            f"Auto-converted {info['converted_count']} gene ID(s) to symbols: "
            f"{', '.join(msg_parts)}."
        )

        if info['failed_count'] > 0:
            conversion_msg += f" ({info['failed_count']} ID(s) could not be converted)"

        # Show some examples
        examples = []
        for orig, sym in (info['ensembl_converted'] + info['entrez_converted'])[:3]:
            examples.append(f"{orig} -> {sym}")
        if examples:
            conversion_msg += f"\nExamples: {', '.join(examples)}"

        print(f"\n[CASSIA] {conversion_msg}\n")

    return converted_markers, info


def is_mygene_available() -> bool:
    """Check if mygene package is available for gene ID conversion."""
    return MYGENE_AVAILABLE


def convert_dataframe_gene_ids(
    df,
    gene_column: str,
    species: str = 'human',
    auto_convert: bool = True,
    verbose: bool = True
) -> Tuple['pd.DataFrame', Dict[str, any]]:
    """
    Convert Ensembl/Entrez gene IDs to symbols in a DataFrame column.

    This is the main entry point for batch-level gene ID conversion.
    Converts IDs once for the entire DataFrame, avoiding repeated API calls
    for each cluster.

    Args:
        df: pandas DataFrame containing gene markers
        gene_column: Name of the column containing gene markers (comma-separated)
        species: Species for gene lookup ('human' or 'mouse')
        auto_convert: Whether to attempt automatic conversion
        verbose: Whether to print conversion summary

    Returns:
        Tuple of:
            - DataFrame with converted gene IDs
            - Conversion info dict with summary statistics
    """
    import pandas as pd

    info = {
        'converted_count': 0,
        'failed_count': 0,
        'total_unique_genes': 0,
        'ensembl_converted': [],
        'entrez_converted': [],
        'failed_ids': [],
        'conversion_performed': False
    }

    if not auto_convert:
        return df, info

    # Make a copy to avoid modifying the original
    df = df.copy()

    # Collect all unique gene IDs from the column
    all_genes = set()
    for value in df[gene_column].dropna():
        if isinstance(value, str):
            # Split by comma and clean
            genes = [g.strip() for g in value.split(',') if g.strip()]
            all_genes.update(genes)

    if not all_genes:
        return df, info

    info['total_unique_genes'] = len(all_genes)

    # Classify all genes
    classified = classify_markers(list(all_genes))

    # If no Ensembl or Entrez IDs, return as-is
    if not classified['ensembl'] and not classified['entrez']:
        if verbose:
            print(f"[CASSIA] All {len(all_genes)} genes appear to be gene symbols. No conversion needed.")
        return df, info

    # Warn if mygene not available (Ensembl REST fallback still works for Ensembl IDs)
    if not MYGENE_AVAILABLE:
        if classified['entrez']:
            warnings.warn(
                f"Detected {len(classified['entrez'])} Entrez IDs but 'mygene' package is not installed. "
                "Entrez conversion requires mygene. Install with: pip install mygene\n"
                "Ensembl IDs will still be converted via Ensembl REST API.",
                UserWarning
            )
        if not classified['ensembl']:
            return df, info

    info['conversion_performed'] = True

    if verbose:
        print(f"\n[CASSIA] Detected gene IDs requiring conversion:")
        if classified['ensembl']:
            print(f"  - {len(classified['ensembl'])} Ensembl IDs")
        if classified['entrez']:
            print(f"  - {len(classified['entrez'])} Entrez IDs")
        print(f"  - {len(classified['symbol'])} gene symbols (no conversion needed)")
        print("Converting IDs to gene symbols...")

    # Convert Ensembl IDs (handles version stripping, species detection, REST fallback)
    ensembl_mapping, ensembl_failed = convert_ensembl_to_symbols(
        classified['ensembl'], species
    )

    # Convert Entrez IDs
    entrez_mapping, entrez_failed = convert_entrez_to_symbols(
        classified['entrez'], species
    )

    # Build combined mapping
    id_to_symbol = {}
    id_to_symbol.update(ensembl_mapping)
    id_to_symbol.update(entrez_mapping)

    # Track conversion stats
    for orig, sym in ensembl_mapping.items():
        info['ensembl_converted'].append((orig, sym))
        info['converted_count'] += 1

    for orig, sym in entrez_mapping.items():
        info['entrez_converted'].append((orig, sym))
        info['converted_count'] += 1

    info['failed_ids'] = ensembl_failed + entrez_failed
    info['failed_count'] = len(info['failed_ids'])

    # Apply mapping to the DataFrame column
    def convert_gene_string(gene_string):
        if not isinstance(gene_string, str):
            return gene_string
        genes = [g.strip() for g in gene_string.split(',') if g.strip()]
        converted = []
        for gene in genes:
            if gene in id_to_symbol:
                converted.append(id_to_symbol[gene])
            else:
                converted.append(gene)
        return ', '.join(converted)

    df[gene_column] = df[gene_column].apply(convert_gene_string)

    # Print summary
    if verbose and info['converted_count'] > 0:
        msg_parts = []
        if info['ensembl_converted']:
            msg_parts.append(f"{len(info['ensembl_converted'])} Ensembl IDs")
        if info['entrez_converted']:
            msg_parts.append(f"{len(info['entrez_converted'])} Entrez IDs")

        print(f"[CASSIA] Successfully converted {info['converted_count']} gene ID(s): {', '.join(msg_parts)}")

        if info['failed_count'] > 0:
            print(f"[CASSIA] Warning: {info['failed_count']} ID(s) could not be converted: "
                  f"{', '.join(info['failed_ids'][:5])}{'...' if len(info['failed_ids']) > 5 else ''}")

        # Show examples
        examples = []
        for orig, sym in (info['ensembl_converted'] + info['entrez_converted'])[:3]:
            examples.append(f"{orig} -> {sym}")
        if examples:
            print(f"[CASSIA] Examples: {', '.join(examples)}\n")

    return df, info
