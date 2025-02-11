# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Sat Dec  7 01:23:04 2024
"""

import os
import requests
import zipfile
import pandas as pd
import matplotlib.pyplot as plt

# Step 1: Download and extract SSA data
def download_and_extract_ssa_data(url, folder="ssa_data"):
    os.makedirs(folder, exist_ok=True)
    zip_path = os.path.join(folder, "names.zip")
    
    # Check if the data has already been downloaded
    if os.path.exists(zip_path):
        print("Data already downloaded.")
    else:
        # Download the zip file
        print("Downloading data...")
        response = requests.get(url)
        if response.status_code == 200:
            with open(zip_path, "wb") as file:
                file.write(response.content)
            print("Data downloaded successfully.")
        else:
            print("Failed to download data.")
            return
    
    # Check if the data has already been extracted
    if len(os.listdir(folder)) > 1:  # If folder contains more than just the zip file
        print("Data already extracted.")
    else:
        # Extract the zip file
        print("Extracting data...")
        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            zip_ref.extractall(folder)
        print("Data extracted successfully.")


# Step 2: Load and combine data from multiple years
def load_ssa_data(start_year, end_year, folder="ssa_data"):
    dataframes = []
    for year in range(start_year, end_year + 1):
        file_path = os.path.join(folder, f"yob{year}.txt")
        if os.path.exists(file_path):
            df = pd.read_csv(file_path, header=None, names=["Name", "Gender", "Count"])
            df["Year"] = year
            dataframes.append(df)
    return pd.concat(dataframes, ignore_index=True) if dataframes else None

# Step 3: Analyze name trends
def analyze_name_trends(data, name):
    filtered_data = data[data["Name"].str.lower() == name.lower()]
    if filtered_data.empty:
        print(f"No data found for the name '{name}'")
        return None, None
    
    # Group and aggregate data
    yearly_counts = filtered_data.groupby(["Year", "Gender"])["Count"].sum().unstack(fill_value=0)
    
    # Ensure missing columns are added with zero values
    if "M" not in yearly_counts.columns:
        yearly_counts["M"] = 0
    if "F" not in yearly_counts.columns:
        yearly_counts["F"] = 0
    
    yearly_counts["Total"] = yearly_counts.sum(axis=1)
    yearly_counts["Male %"] = (yearly_counts["M"] / yearly_counts["Total"]) * 100
    yearly_counts["Female %"] = (yearly_counts["F"] / yearly_counts["Total"]) * 100
    
    return yearly_counts, yearly_counts[["Male %", "Female %"]].mean()


# Step 4: Visualize name trends
def plot_name_trends(yearly_counts, name):
    if yearly_counts is None or yearly_counts.empty:
        print("No trends to plot.")
        return
    
    plt.figure(figsize=(8, 6))  # Increase figure size
    plt.plot(yearly_counts.index, yearly_counts["M"], label="Male", marker="o")
    plt.plot(yearly_counts.index, yearly_counts["F"], label="Female", marker="x")
    plt.plot(yearly_counts.index, yearly_counts["Total"], label="Total", linestyle="--", color="gray")
    
    plt.title(f"Popularity of the Name '{name}' Over Time")
    plt.xlabel("Year")
    plt.ylabel("Number of Babies")
    plt.xticks(ticks=yearly_counts.index[::2], rotation=45)  # Rotate and reduce tick frequency
    plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.tight_layout()  # Adjust layout to fit everything
    plt.show()


from difflib import get_close_matches

# Step 3: Modify the name search logic
def find_closest_name(data, target_name):
    unique_names = data["Name"].unique()
    closest_matches = get_close_matches(target_name, unique_names, n=1, cutoff=0.6)  # Top 3 matches
    if closest_matches:
        print(f"Closest matches to '{target_name}': {closest_matches}")
        return closest_matches
    else:
        print(f"No close matches found for '{target_name}'.")
        return []

# Main Program
if __name__ == "__main__":
    # Define parameters
    start_year = 1980
    end_year = 2022
    
    # Step 1: Input name
    target_name = 'ew2'
    
    # Step 2: Download and extract SSA data
    ssa_url = "https://www.ssa.gov/OACT/babynames/names.zip"
    download_and_extract_ssa_data(ssa_url)
    
    # Step 3: Load data
    ssa_data = load_ssa_data(start_year, end_year)
    if ssa_data is None:
        print("No data loaded.")
    else:
        # Step 4: Find closest names
        closest_names = find_closest_name(ssa_data, target_name)
        
        if closest_names:
            for name in closest_names:
                print(f"\nAnalyzing trends for '{name}'")
                # Step 5: Analyze trends for each close match
                yearly_counts, gender_percentages = analyze_name_trends(ssa_data, name)
                if yearly_counts is not None:
                    print(yearly_counts)
                    print(f"Average Male %: {gender_percentages['Male %']:.2f}%")
                    print(f"Average Female %: {gender_percentages['Female %']:.2f}%")
                    
                    # Step 6: Plot trends
                    plot_name_trends(yearly_counts, name)
