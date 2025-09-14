#Let's import the packages that we'll need for our script
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import os

#Let's also import our own module
import mcb201b_modules.mcb201b_swarmplot_annotated as mcb201b

#Let's move into our small_mol_mtt_data directory
os.chdir('./data/small_mol_mtt_data')

#Let's make some directories that we'll output to
try:
    os.mkdir('analysis_results')
except FileExistsError:
    pass

#Let's pull in our file names
file_names = [name for name in os.listdir() if '.csv' in name]
file_names.sort()

#Let's set up a few lists to hold outputs
compound_hits = []
change_hits = []

#We'll now set up our for-loop:
for i in file_names:
    #First thing that we'll do is set up our base name
    base_name = i.replace('_MTT.csv', '')

    #We can have Python talk to us :)
    print(f'Running analysis for {base_name} MTT dataset')

    #Now we can read in our file
    data = pd.read_csv(i, index_col='Unnamed: 0')

    #Update our DataFrame
    data.rename(index={'A': 'Untreated',
                       'B': f'{base_name}',
                       'C': 'Background'},
                inplace=True,
               )

    #To make the rest of this easier, transpose the DataFrame
    data = data.transpose()

    #Calculate the mean background
    mean_bg = np.mean(data['Background'])

    #Subtract the  mean background from all of our wells
    #Subtractive assignment operation
    data -= mean_bg

    #Find the mean for the control corrected absorbances
    ctrl_mean = np.mean(data['Untreated'])

    #Normalize all the values to the control corrected absorbances
    data_norm = data[['Untreated', f'{base_name}']] / ctrl_mean

    #Let's now run our statistical test
    analysis_result = stats.ttest_ind(data_norm['Untreated'],
                                      data_norm[f'{base_name}'],
                                      equal_var=False,
                                      nan_policy='omit'
                                     )

    #Export our normalized data for reference
    data_norm.to_csv(f'./analysis_results/{base_name}_mtt_processed.csv',
                     index=False
                    )

    #Calculate some descriptive statistics
    #save the output for later
    ctrl_norm_mean = np.mean(data_norm['Untreated'])
    treat_norm_mean = np.mean(data_norm[f'{base_name}'])

    #Difference in the means
    quantified = ctrl_norm_mean - treat_norm_mean
    
    #Let's create an output summary log
    with open(f'./analysis_results/{base_name}_results.txt', 'w') as analysis_file:
        #Write to our file - what our means are
        analysis_file.write(f'Mean of the untreated samples = {ctrl_norm_mean:.2f}\nMean of the treated samples = {treat_norm_mean:.2f}\n')

        analysis_file.write(f'p-value={analysis_result.pvalue:.4f}')

    if analysis_result.pvalue < 0.05:
        compound_hits.append(base_name)
        change_hits.append(quantified)
    else:
        pass

    #Talk to us to see what's going on
    print(f'Generating plot for {base_name}')

    #Let's pull in our module to use for plotting
    mcb201b.swarmplot_annotation(dataset=data_norm,
                                 col_list=data_norm.columns,
                                 conclusion=analysis_result,
                                 y_axislabel='Relative corrected absorbance',
                                 file_name=f'{base_name}_MTT_swarmplot'
                                )

#Talk a little bit more...
print('Generating our list of candidates...')

candidate_hits = pd.DataFrame({'compound': compound_hits,
                               'quantification': change_hits}
                             )

candidate_hits.to_csv('./analysis_results/candidate_hits.csv',
                      index=False
                     )


        














