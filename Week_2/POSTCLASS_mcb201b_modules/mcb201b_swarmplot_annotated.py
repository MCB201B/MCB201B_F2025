#Import needed packages for our module
import numpy 
import pandas
import matplotlib.pyplot
import seaborn
import os

#Define our function for our module
def swarmplot_annotation(dataset, col_list, conclusion, y_axislabel, file_name):
    
    fig, ax = matplotlib.pyplot.subplots()
    dot_colors = ['#777777', '#E64B35']
    
    seaborn.swarmplot(data=dataset[col_list],
                  s=6,
                  palette=dot_colors,
                  zorder=0
                 )
    
    seaborn.barplot(data=dataset[col_list],
                estimator=numpy.mean,
                alpha=0,
                errorbar='se',
                capsize=0.3,
                err_kws={'linewidth': 1,
                         'color': 'k'}
               )
    
    seaborn.boxplot(data=dataset[col_list],
                showmeans=True,
                meanline=True,
                width=0.5,
                meanprops={'linewidth': 1,
                           'color': 'k',
                           'linestyle': '-'
                          },
                medianprops={'visible': False},
                whiskerprops={'visible': False},
                showfliers=False,
                showbox=False,
                showcaps=False
               )
    
    #####################################
    x1, x2 = 0, 1
    y_max = numpy.max(dataset[col_list].max())
    adjustment = y_max*1.1
    h = y_max*0.05
    
    matplotlib.pyplot.plot([x1, x1, x2, x2], [adjustment, adjustment + h, adjustment + h, adjustment], color='k')
    #####################################
    
    #####################################
    if conclusion.pvalue < 0.0005:
        pvalue_star = '***'
    elif conclusion.pvalue < 0.005:
        pvalue_star = '**'
    elif conclusion.pvalue < 0.05:
        pvalue_star = '*'
    else:
        pvalue_star = 'ns'
    
    matplotlib.pyplot.text((x1+x2)/2 , adjustment+(h*1.25), pvalue_star, ha='center', va='bottom', color='k', size=12)
    #####################################
    
    ax.set_ylabel(y_axislabel,
                  size=14
                 )
    
    ax.set_xticks([0,1],
                  col_list,  #Make a little bit more generalizable
                  rotation=35,
                  ha='right',
                  size=12
                 )
    matplotlib.pyplot.yticks(size=12)
    
    fig.set_size_inches(3, 4)
    fig.set_dpi(300)
    
    seaborn.despine()
    matplotlib.pyplot.show()

    #Let's make a directory for our plot
    try:
        os.mkdir('plots')
    except FileExistsError:
        pass

    #We'll save our plot to our 'plots' directory
    fig.savefig(f'./plots/{file_name}.pdf', bbox_inches='tight')

    return None