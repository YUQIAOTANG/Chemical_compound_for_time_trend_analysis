import pandas as pd

# read the data
file1_path = "cleandata_artificial.xlsx"
file2_path = "artificialgroup.xlsx"

cleandata_df = pd.read_excel(file1_path)
artificialgroup_df = pd.read_excel(file2_path)

# remove the space of the data
cleandata_df.columns = cleandata_df.columns.str.strip()
artificialgroup_df.columns = artificialgroup_df.columns.str.strip()

# convert cleandata_df
cleandata_long = cleandata_df.melt(id_vars=["name"], var_name="Condition", value_name="Concentration")

# convert artificialgroup_df
group_long = artificialgroup_df.melt(var_name="Group", value_name="Condition")

# normarized the format of the dataset
def clean_string(s):
    if isinstance(s, str):
        return s.strip().replace("\n", "").replace("\t", "").lower()
    return s

cleandata_long["Condition"] = cleandata_long["Condition"].apply(clean_string)
group_long["Condition"] = group_long["Condition"].apply(clean_string)

# confirm the format are the same
cleandata_long["Condition"] = cleandata_long["Condition"].astype(str)
group_long["Condition"] = group_long["Condition"].astype(str)

# merge the data
merged_data = cleandata_long.merge(group_long, on="Condition", how="left")

# find the unmatched Condition
unmatched_conditions = merged_data[merged_data["Group"].isna()]["Condition"].unique()
if len(unmatched_conditions) > 0:
    print("未匹配的 Condition：")
    print(unmatched_conditions)

# save the new datasets
merged_data.to_excel("grouped_data_fixed.xlsx", index=False)
print("data classify has been done, save as grouped_data_fixed.xlsx")
