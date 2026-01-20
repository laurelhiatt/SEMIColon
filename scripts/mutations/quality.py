import json
import os
import glob
import csv
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Path to your directory with the JSON files
folder_path = '/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/reports'

# Prepare containers
laurel_q20 = []
laurel_q30 = []
other_q20 = []
other_q30 = []

hiatt_donors = ["AS", "AC", "DC", "DE", "CE", "RE", "TR", "SI", "TC", "RTM", "Laurel", "GIB"]


def filename_contains_any_string(filename, string_list):
    """
    Checks if a filename contains any of the strings in a given list.

    Args:
        filename (str): The name of the file to check.
        string_list (list): A list of strings to search for within the filename.

    Returns:
        bool: True if the filename contains any string from the list, False otherwise.
    """
    for search_string in string_list:
        if search_string in filename:
            return True
    return False

# Prepare a list for CSV export
rows = []

# Go through each fastp report
for filepath in glob.glob(os.path.join(folder_path, '*-fastp-report.json')):
    with open(filepath, 'r') as f:
        data = json.load(f)

        # Extract rates
        q20_rate = data['summary']['before_filtering']['q20_rate']
        q30_rate = data['summary']['before_filtering']['q30_rate']

        # Extract sample name (removing '-fastp-report.json' part)
        filename = os.path.basename(filepath)
        sample_name = filename.replace('-fastp-report.json', '')

        # Decide group
        if filename_contains_any_string(filename, hiatt_donors):
            group = 'Laurel'
            laurel_q20.append(q20_rate)
            laurel_q30.append(q30_rate)
        else:
            group = 'Other'
            other_q20.append(q20_rate)
            other_q30.append(q30_rate)

        # Append a row for the CSV
        rows.append({
            'sample_name': sample_name,
            'group': group,
            'q20_rate': q20_rate,
            'q30_rate': q30_rate
        })

# Function to safely compute averages
def average(lst):
    return sum(lst) / len(lst) if lst else float('nan')

# Print results
print(f"Laurel samples: {len(laurel_q20)}")
print(f"  Avg Q20 rate: {average(laurel_q20):.4f}")
print(f"  Avg Q30 rate: {average(laurel_q30):.4f}")

print(f"Other samples: {len(other_q20)}")
print(f"  Avg Q20 rate: {average(other_q20):.4f}")
print(f"  Avg Q30 rate: {average(other_q30):.4f}")

# Save to CSV
output_csv = os.path.join(folder_path, 'fastp_rates_summary.csv')
with open(output_csv, 'w', newline='') as csvfile:
    fieldnames = ['sample_name', 'group', 'q20_rate', 'q30_rate']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    writer.writerows(rows)

print(f"\nSummary CSV written to: {output_csv}")

# Load the CSV we just created
df = pd.read_csv(output_csv)

# Set plot style
sns.set(style="whitegrid")

# Melt the DataFrame to long format so we can plot Q20 and Q30 together
df_melted = df.melt(id_vars=['sample_name', 'group'], value_vars=['q20_rate', 'q30_rate'],
                    var_name='Quality Metric', value_name='Rate')

# Create the plot
plt.figure(figsize=(10,6))
sns.boxplot(data=df_melted, x='Quality Metric', y='Rate', hue='group')

# Tweak labels and title
plt.title('Q20 and Q30 Rates by Sample Group')
plt.xlabel('Quality Metric')
plt.ylabel('Rate')
plt.ylim(0.7, 1.0)  # Optional: zoom in to typical quality score range
plt.legend(title='Sample Group')

# Show it
plt.tight_layout()

# Save the plot
plot_path = os.path.join(folder_path, 'fastp_rates_boxplot.png')
plt.savefig(plot_path, dpi=300)
print(f"Plot saved to {plot_path}")

plt.show()
