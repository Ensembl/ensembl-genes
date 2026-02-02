"""
BUSCO Results Comparison Tool
Compares two BUSCO outputs and identifies which genes are missed in one annotation
"""
import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

def parse_busco_full_table(busco_file):
    """Parse BUSCO full_table.tsv file"""
    buscos = {}

    with open(busco_file, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue

            busco_id = parts[0]
            status = parts[1]  # Complete, Duplicated, Fragmented, Missing

            # Store additional info if available
            info = {
                'status': status,
                'sequence': parts[2] if len(parts) > 2 else None,
                'score': parts[3] if len(parts) > 3 else None,
                'length': parts[4] if len(parts) > 4 else None
            }

            buscos[busco_id] = info

    return buscos

def categorize_buscos(buscos):
    """Categorize BUSCOs by status"""
    categories = {
        'Complete': set(),
        'Duplicated': set(),
        'Fragmented': set(),
        'Missing': set()
    }

    for busco_id, info in buscos.items():
        status = info['status']
        # Handle both "Complete" and "Duplicated" status
        if status == 'Complete':
            categories['Complete'].add(busco_id)
        elif status == 'Duplicated':
            categories['Duplicated'].add(busco_id)
            categories['Complete'].add(busco_id)  # Duplicated are also complete
        elif status == 'Fragmented':
            categories['Fragmented'].add(busco_id)
        elif status == 'Missing':
            categories['Missing'].add(busco_id)

    return categories

def compare_busco_results(file1, file2, name1="Funannotate", name2="EnsemblAnno"):
    """Compare two BUSCO results and identify differences"""

    print(f"Parsing {name1} BUSCO results...")
    buscos1 = parse_busco_full_table(file1)

    print(f"Parsing {name2} BUSCO results...")
    buscos2 = parse_busco_full_table(file2)

    # Categorize BUSCOs
    cat1 = categorize_buscos(buscos1)
    cat2 = categorize_buscos(buscos2)

    # Calculate actual complete counts (including duplicated)
    complete1 = len([b for b in buscos1.values() if b['status'] in ['Complete', 'Duplicated']])
    complete2 = len([b for b in buscos2.values() if b['status'] in ['Complete', 'Duplicated']])
    single1 = len([b for b in buscos1.values() if b['status'] == 'Complete'])
    single2 = len([b for b in buscos2.values() if b['status'] == 'Complete'])
    dup1 = len([b for b in buscos1.values() if b['status'] == 'Duplicated'])
    dup2 = len([b for b in buscos2.values() if b['status'] == 'Duplicated'])
    frag1 = len(cat1['Fragmented'])
    frag2 = len(cat2['Fragmented'])
    miss1 = len(cat1['Missing'])
    miss2 = len(cat2['Missing'])

    # Print summary statistics
    print("\n" + "="*80)
    print("BUSCO SUMMARY STATISTICS")
    print("="*80)

    total1 = len(buscos1)
    total2 = len(buscos2)

    print(f"\n{'Category':<25} {name1:>20} {name2:>20} {'Difference':>12}")
    print("-"*80)

    print(f"{'Complete (C)':<25} {complete1:>10} ({complete1/total1*100:>5.1f}%) {complete2:>10} ({complete2/total2*100:>5.1f}%) {complete1-complete2:>12}")
    print(f"  {'Single copy (S)':<23} {single1:>10} ({single1/total1*100:>5.1f}%) {single2:>10} ({single2/total2*100:>5.1f}%) {single1-single2:>12}")
    print(f"  {'Duplicated (D)':<23} {dup1:>10} ({dup1/total1*100:>5.1f}%) {dup2:>10} ({dup2/total2*100:>5.1f}%) {dup1-dup2:>12}")
    print(f"{'Fragmented (F)':<25} {frag1:>10} ({frag1/total1*100:>5.1f}%) {frag2:>10} ({frag2/total2*100:>5.1f}%) {frag1-frag2:>12}")
    print(f"{'Missing (M)':<25} {miss1:>10} ({miss1/total1*100:>5.1f}%) {miss2:>10} ({miss2/total2*100:>5.1f}%) {miss1-miss2:>12}")
    print(f"{'Total BUSCOs (n)':<25} {total1:>20} {total2:>20} {total1-total2:>12}")

    # Verify percentages match what user reported
    print(f"\n{name1}: C:{complete1/total1*100:.1f}%[S:{single1/total1*100:.1f}%,D:{dup1/total1*100:.1f}%],F:{frag1/total1*100:.1f}%,M:{miss1/total1*100:.1f}%,n:{total1}")
    print(f"{name2}: C:{complete2/total2*100:.1f}%[S:{single2/total2*100:.1f}%,D:{dup2/total2*100:.1f}%],F:{frag2/total2*100:.1f}%,M:{miss2/total2*100:.1f}%,n:{total2}")

    # Analyze what Anno is missing
    print("\n" + "="*80)
    print(f"DETAILED ANALYSIS: What {name2} is missing compared to {name1}")
    print("="*80)

    # BUSCOs that are Complete in name1 but not in name2
    complete_busco_ids_1 = {bid for bid, info in buscos1.items() if info['status'] in ['Complete', 'Duplicated']}
    complete_busco_ids_2 = {bid for bid, info in buscos2.items() if info['status'] in ['Complete', 'Duplicated']}

    complete_in_1_not_2 = complete_busco_ids_1 - complete_busco_ids_2

    # Check what happened to these BUSCOs in name2
    status_breakdown = defaultdict(list)

    for busco_id in complete_in_1_not_2:
        status_in_2 = buscos2.get(busco_id, {}).get('status', 'Not found')
        status_breakdown[status_in_2].append(busco_id)

    print(f"\n{name1} has {complete1} Complete BUSCOs (Single + Duplicated)")
    print(f"{name2} has {complete2} Complete BUSCOs (Single + Duplicated)")
    print(f"\nBUSCOs Complete in {name1} but NOT Complete in {name2}: {len(complete_in_1_not_2)}")

    if len(complete_in_1_not_2) > 0:
        print(f"\nStatus of these {len(complete_in_1_not_2)} BUSCOs in {name2}:")
        print("-"*50)
        for status, busco_list in sorted(status_breakdown.items(), key=lambda x: -len(x[1])):
            print(f"  {status:<20}: {len(busco_list):>5} BUSCOs ({len(busco_list)/len(complete_in_1_not_2)*100:.1f}%)")

            # Print first few examples
            if len(busco_list) <= 10:
                for busco in busco_list:
                    print(f"    - {busco}")
            else:
                for busco in busco_list[:5]:
                    print(f"    - {busco}")
                print(f"    ... and {len(busco_list) - 5} more")

    # BUSCOs that are Missing in name2 but present in name1
    print("\n" + "="*80)
    print(f"BUSCOs that are Missing in {name2} but found in {name1}:")
    print("="*80)

    missing_in_2_ids = {bid for bid, info in buscos2.items() if info['status'] == 'Missing'}
    found_in_1 = set()

    status_in_1 = defaultdict(list)
    for busco_id in missing_in_2_ids:
        if busco_id in buscos1 and buscos1[busco_id]['status'] != 'Missing':
            found_in_1.add(busco_id)
            status_in_1[buscos1[busco_id]['status']].append(busco_id)

    print(f"\n{name2} is Missing {miss2} BUSCOs")
    print(f"Of these, {len(found_in_1)} are found (not missing) in {name1}")

    if len(found_in_1) > 0:
        print(f"\nStatus of these BUSCOs in {name1}:")
        print("-"*50)
        for status, busco_list in sorted(status_in_1.items(), key=lambda x: -len(x[1])):
            print(f"  {status:<20}: {len(busco_list):>5} BUSCOs ({len(busco_list)/len(found_in_1)*100:.1f}%)")

            # Print examples
            if len(busco_list) <= 10:
                for busco in busco_list:
                    print(f"    - {busco}")
            else:
                for busco in busco_list[:5]:
                    print(f"    - {busco}")
                print(f"    ... and {len(busco_list) - 5} more")

    # Fragmented analysis
    print("\n" + "="*80)
    print(f"Fragmented BUSCOs analysis:")
    print("="*80)

    frag_in_2_ids = {bid for bid, info in buscos2.items() if info['status'] == 'Fragmented'}
    complete_in_1 = set()

    for busco_id in frag_in_2_ids:
        if busco_id in buscos1 and buscos1[busco_id]['status'] in ['Complete', 'Duplicated']:
            complete_in_1.add(busco_id)

    print(f"\n{name2} has {frag2} Fragmented BUSCOs")
    print(f"Of these, {len(complete_in_1)} are Complete in {name1}")

    if len(complete_in_1) > 0:
        print(f"\nThese {len(complete_in_1)} BUSCOs are Fragmented in {name2} but Complete in {name1}:")
        if len(complete_in_1) <= 20:
            for busco in sorted(complete_in_1):
                print(f"  - {busco}")
        else:
            for busco in sorted(list(complete_in_1)[:10]):
                print(f"  - {busco}")
            print(f"  ... and {len(complete_in_1) - 10} more")

    # Create visualizations
    create_comparison_plots(complete1, single1, dup1, frag1, miss1,
                           complete2, single2, dup2, frag2, miss2,
                           name1, name2, complete_in_1_not_2,
                           status_breakdown, found_in_1, status_in_1)

    # Save detailed results to file
    save_detailed_results(complete_in_1_not_2, status_breakdown, found_in_1,
                         status_in_1, complete_in_1, name1, name2)

    return buscos1, buscos2, complete_in_1_not_2, missing_in_2_ids

def create_comparison_plots(complete1, single1, dup1, frag1, miss1,
                           complete2, single2, dup2, frag2, miss2,
                           name1, name2, complete_in_1_not_2,
                           status_breakdown, found_in_1, status_in_1):
    """Create visualization plots for BUSCO comparison"""

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Overall BUSCO status comparison (showing single, duplicated separately)
    categories = ['Single\nCopy', 'Duplicated', 'Fragmented', 'Missing']
    counts1 = [single1, dup1, frag1, miss1]
    counts2 = [single2, dup2, frag2, miss2]

    x = range(len(categories))
    width = 0.35

    axes[0, 0].bar([i - width/2 for i in x], counts1, width, label=name1,
                   color='steelblue', alpha=0.8)
    axes[0, 0].bar([i + width/2 for i in x], counts2, width, label=name2,
                   color='coral', alpha=0.8)
    axes[0, 0].set_ylabel('Number of BUSCOs')
    axes[0, 0].set_title('BUSCO Status Comparison')
    axes[0, 0].set_xticks(x)
    axes[0, 0].set_xticklabels(categories, rotation=0, ha='center')
    axes[0, 0].legend()
    axes[0, 0].grid(axis='y', alpha=0.3)

    # Plot 2: Percentage stacked bar chart
    total1 = complete1 + frag1 + miss1
    total2 = complete2 + frag2 + miss2

    pct1_complete = complete1/total1*100
    pct1_frag = frag1/total1*100
    pct1_miss = miss1/total1*100

    pct2_complete = complete2/total2*100
    pct2_frag = frag2/total2*100
    pct2_miss = miss2/total2*100

    annotations = [name1, name2]
    complete = [pct1_complete, pct2_complete]
    fragmented = [pct1_frag, pct2_frag]
    missing = [pct1_miss, pct2_miss]

    axes[0, 1].barh(annotations, complete, label='Complete', color='#2ecc71')
    axes[0, 1].barh(annotations, fragmented, left=complete, label='Fragmented', color='#f39c12')
    axes[0, 1].barh(annotations, missing,
                    left=[complete[i]+fragmented[i] for i in range(2)],
                    label='Missing', color='#e74c3c')

    # Add percentage labels
    for i, ann in enumerate(annotations):
        axes[0, 1].text(complete[i]/2, i, f'{complete[i]:.1f}%',
                       ha='center', va='center', fontweight='bold', color='white')
        axes[0, 1].text(complete[i] + fragmented[i]/2, i, f'{fragmented[i]:.1f}%',
                       ha='center', va='center', fontweight='bold')
        axes[0, 1].text(complete[i] + fragmented[i] + missing[i]/2, i, f'{missing[i]:.1f}%',
                       ha='center', va='center', fontweight='bold', color='white')

    axes[0, 1].set_xlabel('Percentage (%)')
    axes[0, 1].set_title('BUSCO Status Distribution (%)')
    axes[0, 1].legend(loc='lower right')
    axes[0, 1].set_xlim(0, 100)

    # Plot 3: What happened to Complete BUSCOs from name1 in name2
    if len(complete_in_1_not_2) > 0:
        statuses = list(status_breakdown.keys())
        counts = [len(status_breakdown[s]) for s in statuses]
        colors = {'Missing': '#e74c3c', 'Fragmented': '#f39c12',
                 'Duplicated': '#3498db', 'Not found': '#95a5a6'}
        bar_colors = [colors.get(s, '#7f8c8d') for s in statuses]

        axes[1, 0].bar(range(len(statuses)), counts, color=bar_colors, alpha=0.8)
        axes[1, 0].set_ylabel('Number of BUSCOs')
        axes[1, 0].set_title(f'Status in {name2} of BUSCOs\nComplete in {name1} but not in {name2}\n(n={len(complete_in_1_not_2)})')
        axes[1, 0].set_xticks(range(len(statuses)))
        axes[1, 0].set_xticklabels(statuses, rotation=15, ha='right')
        axes[1, 0].grid(axis='y', alpha=0.3)

        # Add count labels on bars
        for i, count in enumerate(counts):
            axes[1, 0].text(i, count, str(count), ha='center', va='bottom')
    else:
        axes[1, 0].text(0.5, 0.5, 'No differences found',
                       ha='center', va='center', fontsize=12)
        axes[1, 0].set_title(f'Complete BUSCOs: No differences')

    # Plot 4: Missing in name2 but found in name1
    if len(found_in_1) > 0:
        statuses = list(status_in_1.keys())
        counts = [len(status_in_1[s]) for s in statuses]
        colors = {'Complete': '#2ecc71', 'Fragmented': '#f39c12',
                 'Duplicated': '#3498db'}
        bar_colors = [colors.get(s, '#7f8c8d') for s in statuses]

        axes[1, 1].bar(range(len(statuses)), counts, color=bar_colors, alpha=0.8)
        axes[1, 1].set_ylabel('Number of BUSCOs')
        axes[1, 1].set_title(f'Status in {name1} of BUSCOs\nMissing in {name2}\n(n={len(found_in_1)})')
        axes[1, 1].set_xticks(range(len(statuses)))
        axes[1, 1].set_xticklabels(statuses, rotation=15, ha='right')
        axes[1, 1].grid(axis='y', alpha=0.3)

        # Add count labels on bars
        for i, count in enumerate(counts):
            axes[1, 1].text(i, count, str(count), ha='center', va='bottom')
    else:
        axes[1, 1].text(0.5, 0.5, 'All Missing BUSCOs\nare also missing in both',
                       ha='center', va='center', fontsize=12)
        axes[1, 1].set_title(f'Missing BUSCOs: Consistent')

    plt.tight_layout()
    plt.savefig('busco_comparison.png', dpi=300, bbox_inches='tight')
    print("\n✓ Plots saved to 'busco_comparison.png'")

def save_detailed_results(complete_in_1_not_2, status_breakdown, found_in_1,
                         status_in_1, complete_in_1, name1, name2):
    """Save detailed results to text files"""

    # Save list of BUSCOs complete in name1 but not in name2
    with open(f'buscos_complete_in_{name1}_not_{name2}.txt', 'w') as f:
        f.write(f"# BUSCOs that are Complete in {name1} but NOT Complete in {name2}\n")
        f.write(f"# Total: {len(complete_in_1_not_2)}\n\n")

        for status, busco_list in sorted(status_breakdown.items(), key=lambda x: -len(x[1])):
            f.write(f"\n## Status in {name2}: {status} ({len(busco_list)} BUSCOs)\n")
            for busco in sorted(busco_list):
                f.write(f"{busco}\n")

    print(f"✓ Saved detailed list to 'buscos_complete_in_{name1}_not_{name2}.txt'")

    # Save list of BUSCOs missing in name2 but found in name1
    if len(found_in_1) > 0:
        with open(f'buscos_missing_in_{name2}_found_in_{name1}.txt', 'w') as f:
            f.write(f"# BUSCOs that are Missing in {name2} but found in {name1}\n")
            f.write(f"# Total: {len(found_in_1)}\n\n")

            for status, busco_list in sorted(status_in_1.items(), key=lambda x: -len(x[1])):
                f.write(f"\n## Status in {name1}: {status} ({len(busco_list)} BUSCOs)\n")
                for busco in sorted(busco_list):
                    f.write(f"{busco}\n")

        print(f"✓ Saved detailed list to 'buscos_missing_in_{name2}_found_in_{name1}.txt'")

    # Save fragmented analysis
    if len(complete_in_1) > 0:
        with open(f'buscos_fragmented_in_{name2}_complete_in_{name1}.txt', 'w') as f:
            f.write(f"# BUSCOs that are Fragmented in {name2} but Complete in {name1}\n")
            f.write(f"# Total: {len(complete_in_1)}\n\n")
            for busco in sorted(complete_in_1):
                f.write(f"{busco}\n")

        print(f"✓ Saved detailed list to 'buscos_fragmented_in_{name2}_complete_in_{name1}.txt'")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="Compare two genome BUSCO results for different tools."
	)

	parser.add_argument(
		"--file1",
		required=True,
		help="Path to the 'full_table.tsv' files from BUSCO output directory tool1"
	)
	parser.add_argument(
		"--file2",
		required=True,
		help="Path to the 'full_table.tsv' files from BUSCO output directory tool2"
	)
	parser.add_argument(
		"--name1",
		default="Funannotate",
		help="Name of first annotation tool"
	)
	parser.add_argument(
		"--name2",
		default="EnsemblAnno",
		help="Name of second annotation tool"
	)

	args = parser.parse_args()

	compare_busco_results(args.file1, args.file2, args.name1, args.name2)