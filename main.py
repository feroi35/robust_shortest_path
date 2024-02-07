import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math


def merge_dataframe():
    """
    select two dataframe of simulations, and merge them into a single dataframe, chosing the best objective and lower bound
    for each instance
    """
    df = pd.read_csv('results/dualized_results.csv')
    df2 = pd.read_csv('results/dualized_results_bis.csv')

    df_merge = pd.DataFrame(columns = df.columns)

    for index, row in df.iterrows():
        if row['objective'] != row['lower_bound']:                
            row_retested = df2[df2['instance'].str.contains(row['instance'][-22:])]
            row_retested.at[row_retested.index[0], 'instance'] = row['instance']
            if math.isnan(row['objective']):
                df_merge = pd.concat([df_merge, row_retested], ignore_index=True)
                continue
            if row_retested['objective'].values[0] - row_retested['lower_bound'].values[0] < row['objective'] - row['lower_bound']:
                df_merge = pd.concat([df_merge, row_retested], ignore_index=True)
            elif row_retested['objective'].values[0] - row_retested['lower_bound'].values[0] == row['objective'] - row['lower_bound'] and row_retested['time'].values[0] < row['time']:
                df_merge = pd.concat([df_merge, row_retested], ignore_index=True)
            else:
                df_merge = pd.concat([df_merge, row.to_frame().T], ignore_index=True)
        else:
            df_merge = pd.concat([df_merge, row.to_frame().T], ignore_index=True)
    
    df_merge.to_csv('results/dualized_results_merged.csv', index=False)


def identify_unsolved_instances():
    """
    find the instances where an error happened and ,,,,, was returned
    """
    df = pd.read_csv('results/branch_and_cut_results.csv')
    count = 0
    for index, row in df.iterrows():
        if math.isnan(row['objective']):
            print(row['instance'])
            count += 1
    print(count)
    

def main():
    df = pd.read_csv('test.csv')
    print(df.head())

    df['bestObjective'] = df.groupby('instance')['objective'].transform('min')
    df['bestLowerBound'] = df.groupby('instance')['lowerBound'].transform('max')

    df['gap'] = 100*(df['objective'] - df['bestLowerBound']) / df['bestLowerBound']
    df['gap'] = df['gap'].apply(lambda x: min(100, x))

    # df['PR_inf'] = 100*(df['bestLowerBound'] - df['bestObjective']) / df['bestLowerBound']
    # df['PR_inf'] = df['PR_inf'].apply(lambda x: min(100, x))

    # df['PR_sup'] = 100*(df['bestUpperBound'] - df['bestObjective']) / df['bestUpperBound']


    print(df)


    ######################
    methods = ["a", "b", "static"]
    columns = ["instance"]
    for method in methods:
        if method!='static':
            columns.append(f'objective_{method}')
            columns.append(f'gap_{method}')
        else:
            columns.append('lowerBound_static')
            columns.append('bestLowerBound')
            columns.append('bestUpperBound')

    result_df = pd.DataFrame(columns=columns)

    first_method = True
    for instance in df['instance'].unique():
        instance_df = pd.DataFrame({'instance': [instance]})

        for method in methods:
            if method!='static':
                method_df = df[(df['instance'] == instance) & (df['method'] == method)].reset_index(drop=True)
                instance_df[f'objective_{method}'] = method_df['objective'][0]
                instance_df[f'gap_{method}'] = method_df['gap'][0]
            else:
                static_method_df = df[(df['instance'] == instance) & (df['method'] == 'static')].reset_index(drop=True)
                instance_df['lowerBound_static'] = static_method_df['lowerBound'][0]
                instance_df['bestLowerBound'] = df[df['instance'] == instance]['bestLowerBound'].values[0]
                instance_df['bestObjective'] = df[df['instance'] == instance]['bestObjective'].values[0]

        if first_method:
            result_df = instance_df
            first_method = False
        else:
            result_df = pd.concat([result_df, instance_df], ignore_index=True)

    print(result_df.head())


if __name__ == '__main__':
    identify_unsolved_instances()
