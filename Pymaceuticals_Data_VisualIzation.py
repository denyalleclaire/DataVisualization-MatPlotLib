#!/usr/bin/env python
# coding: utf-8

# # Pymaceuticals Inc.
# ---
# 
# ### Analysis
# 
# Based on the given data, it is clear that the more mass (weight) the mouse has, the larger the mass fo the tumor was. Additionally, it appears as though the length of time the treatment is administered directly coorelates to the size of the tumor, as the longer the treatment is given, the more effective it is as evidenced by the tumor size decreasing. 
#  

# In[1]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st

# Study data files
mouse_metadata_path = "./data/Mouse_metadata.csv"
study_results_path = "./data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single DataFrame
study_data_combined_df= pd.merge(study_results, mouse_metadata, how= "left", on="Mouse ID" )

# Display the data table for preview
study_data_combined_df


# In[2]:


# Checking the number of mice.
len(study_data_combined_df["Mouse ID"].unique())


# In[3]:


# Our data should be uniquely identified by Mouse ID and Timepoint
# Get the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
duplicate_mouseID=study_data_combined_df[study_data_combined_df.duplicated(subset=["Mouse ID", "Timepoint"])]["Mouse ID"].unique()
duplicate_mouseID


# In[4]:


# Optional: Get all the data for the duplicate mouse ID. 
dumplicated_mouse_data=study_data_combined_df[study_data_combined_df["Mouse ID"]=="g989"]
dumplicated_mouse_data


# In[5]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean_data_study=study_data_combined_df[study_data_combined_df["Mouse ID"].isin(duplicate_mouseID)== False]
clean_data_study


# In[6]:


# Checking the number of mice in the clean DataFrame.
len(clean_data_study["Mouse ID"].unique())


# ## Summary Statistics

# In[7]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
# Assemble the resulting series into a single summary DataFrame.

means=clean_data_study.groupby("Drug Regimen").mean()["Tumor Volume (mm3)"]
medians=clean_data_study.groupby("Drug Regimen").median()["Tumor Volume (mm3)"]
variances=clean_data_study.groupby("Drug Regimen").var()["Tumor Volume (mm3)"]
standard_dev=clean_data_study.groupby("Drug Regimen").std()["Tumor Volume (mm3)"]
sems=clean_data_study.groupby("Drug Regimen").sem()["Tumor Volume (mm3)"]


summary_table=pd.DataFrame({"Mean Tumor Volume": means, 
             "Median Tumor Volume": medians, 
             "Tumor Volume Variance": variances, 
             "Tumor Volume Standard Dev.": standard_dev,
             "Tumor Volume Standard Err. ": sems})

summary_table


# In[8]:


# A more advanced method to generate a summary statistics table of mean, median, variance, standard deviation,
# and SEM of the tumor volume for each regimen (only one method is required in the solution)

# Using the aggregation method, produce the same summary statistics in a single line
summary_table_2=clean_data_study.groupby("Drug Regimen").agg({"Tumor Volume (mm3)" :["mean", "median","var", "std","sem" ]})

summary_table_2


# ## Bar and Pie Charts

# In[9]:


# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using Pandas.
counts=clean_data_study["Drug Regimen"].value_counts()
counts.plot(kind="bar")
plt.xlabel("Drug Regimen")
plt.xticks(rotation=90)
plt.ylabel("Number of Mice Tested")

plt.show()


# In[10]:


# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using pyplot.
counts= clean_data_study["Drug Regimen"].value_counts()
plt.bar(counts.index.values, counts.values)
plt.xlabel("Drug Regimen")
plt.xticks(rotation=90)
plt.ylabel("Number of Mice Tested")

plt.show()


# In[11]:


# Generate a pie plot showing the distribution of female versus male mice using Pandas
counts=clean_data_study.Sex.value_counts()
counts.plot(kind="pie", autopct="%1.1f%%")


# In[12]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot
counts=clean_data_study.Sex.value_counts()
plt.pie(counts.values, labels=counts.index.values, autopct=("%1.1f%%"))
plt.show()


# ## Quartiles, Outliers and Boxplots

# In[13]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin

# Start by getting the last (greatest) timepoint for each mouse
max_tumor=clean_data_study.groupby(["Mouse ID"])["Timepoint"].max()
max_tumor= max_tumor.reset_index()
# Merge this group df with the original DataFrame to get the tumor volume at the last timepoint
merged_data=max_tumor.merge(clean_data_study, on=["Mouse ID", "Timepoint"], how="left")
merged_data


# In[14]:


# Put treatments into a list for for loop (and later for plot labels)
treatment_list= ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]


# Create empty list to fill with tumor vol data (for plotting)
tumor_vol_list=[]

# Calculate the IQR and quantitatively determine if there are any potential outliers. 
for drug in treatment_list: 
    
    
    # Locate the rows which contain mice on each drug and get the tumor volumes
    tumor_volume_final=merged_data.loc[merged_data["Drug Regimen"]== drug, "Tumor Volume (mm3)"]
    
    # add subset 
    tumor_vol_list.append(tumor_volume_final)
    
    # Determine outliers using upper and lower bounds
    quartiles= tumor_volume_final.quantile([.25, .5, .75])
    upperq= quartiles[0.75]
    lowerq= quartiles[0.25]
    iqr= upperq-lowerq
    lower_bound= lowerq- (1.5* iqr) 
    upper_bound= upperq+ (1.5* iqr) 
    
    outliers= tumor_volume_final.loc[(tumor_volume_final < lower_bound) | (tumor_volume_final > upper_bound)]
    
    print(f" {drug}'s potential outliers {outliers}")


# In[15]:


# Generate a box plot that shows the distrubution of the tumor volume for each treatment group.
orange_out=dict(markerfacecolor= "red", markersize= 15)
plt.boxplot(tumor_vol_list, labels= treatment_list, flierprops=orange_out)
plt.ylabel("Final Tumor Volume (mm3)")
plt.show()


# ## Line and Scatter Plots

# In[26]:


# Generate a line plot of tumor volume vs. time point for a single mouse treated with Capomulin
Camopulin_Table=clean_data_study[clean_data_study["Drug Regimen"]== "Capomulin"]
mouse_data= Camopulin_Table[Camopulin_Table["Mouse ID"] == "b128"]
plt.plot(mouse_data["Timepoint"], mouse_data["Tumor Volume (mm3)"])
plt.xlabel("Timepoint (days)")
plt.ylabel("Tumor Volume (mm3")
plt.title("Capomulin Treatment Of Mouse B128")
mouse_data


# In[28]:


# Generate a scatter plot of mouse weight vs. the average observed tumor volume for the entire Capomulin regimen
clean_data_study[clean_data_study["Drug Regimen"]=="Camopulin"]
Capomulin_average= Camopulin_Table.groupby(["Mouse ID"]).mean()
plt.scatter(Capomulin_average["Weight (g)"], Capomulin_average["Tumor Volume (mm3)"])
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.show()


# ## Correlation and Regression

# In[43]:


# Calculate the correlation coefficient and a linear regression model 
# for mouse weight and average observed tumor volume for the entire Capomulin regimen

corr=st.pearsonr(Capomulin_average["Weight (g)"],Capomulin_average["Tumor Volume (mm3)"] )
print(f" The correlation between mouse weight and the avergae tumor volume is {round(corr[0],2)}")

model=st.linregress(Capomulin_average["Weight (g)"], Capomulin_average["Tumor Volume (mm3)"])
slope= model[0]
b= model[1]
y_values= Capomulin_average["Weight (g)"] * slope + b 
plt.scatter( Capomulin_average["Weight (g)"], Capomulin_average["Tumor Volume (mm3)"])
plt.plot (Capomulin_average["Weight (g)"], y_values, color="green") 
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")

plt.show()


# In[ ]:




