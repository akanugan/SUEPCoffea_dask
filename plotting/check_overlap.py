import re

# Function to extract runs from the text file
def extract_runs_from_text_file(text_file):
    runs = []
    with open(text_file, 'r') as file:
        for line in file:
            match = re.search(r'(\d+)\.root', line)
            if match:
                runs.append(match.group(1))
    return runs

# Function to read runs from the CSV file
def read_runs_from_csv(csv_file):
    runs = []
    with open(csv_file, 'r') as file:
        next(file)  # Skip header
        for line in file:
            run = line.split(',')[0].strip()
            runs.append(run)
    return runs

# Main function to compare runs from both files
def compare_runs(csv_file, text_file):
    csv_runs = read_runs_from_csv(csv_file)
    text_runs = extract_runs_from_text_file(text_file)
    
    matching_runs = set(csv_runs) & set(text_runs)
    non_matching_csv = set(csv_runs) - matching_runs
    non_matching_text = set(text_runs) - matching_runs
    percentage_matched = (len(matching_runs) / len(csv_runs)) * 100
    
    print("Percentage of matched runs:", percentage_matched)
    print("Matching runs:")
    for run in matching_runs:
        print(run)
    
    print("\nNon-matching runs in CSV file:")
    for run in non_matching_csv:
        print(run)
    
    print("\nNon-matching runs in text file:")
    for run in non_matching_text:
        print(run)

    

# Example usage
if __name__ == "__main__":
    csv_file = "plotting/outliers_3e6_3e7_noHT.csv"
    text_file = "plotting/OMSruns.txt"
    compare_runs(csv_file, text_file)