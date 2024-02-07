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


def make_graphic(big_df):
    sns.set()
    fig = plt.figure(figsize=(10, 5))
    time = np.arange(1, 1200, 1)
    for method in big_df['method'].unique():
        print(method)
        nb_closed_instances = []
        # compter le nombre d'instances closed en moins que le temps
        for t in time:
            nb_closed_instances.append(big_df[(big_df['method'] == method) & (big_df['time'] < t) & (big_df['closed'] == True)].shape[0])
        plt.plot(time, nb_closed_instances, label=method)
    plt.xlabel("Time (s)")
    plt.ylabel("Number of closed instances")
    plt.legend()
    plt.savefig("aaaa.png")
    return


def make_results_tab():
    # methods = ["plans_coupants", "branch_and_cut", "dualized", "heuristic", "static"]
    methods = ["plans_coupants", "dualized", "heuristic", "static"]

    # df = pd.read_csv("test.csv")
    # methods = ["a", "b", "static"]

    df_initialized = False
    for method in methods:
        if not df_initialized:
            df = pd.read_csv(f'results/{method}_results.csv')
            df_initialized = True
        else:
            df = pd.concat([df, pd.read_csv(f'results/{method}_results.csv')], ignore_index=True)

    # On recupére la meilleure solution obtenue pour chaque instance à l'aide des différentes méthodes, pour calculer un sup du Poids de la Robustesse
    df['bestObjective'] = df.groupby('instance')['objective'].transform('min')
    # On recupére la meilleure borne inf obtenue pour chaque instance à l'aide des différentes méthodes (souvent dualized)
    df['bestlower_bound'] = df.groupby('instance')['lower_bound'].transform('max')

    df['gap'] = 100*(df['objective'] - df['bestlower_bound']) / df['bestlower_bound']
    df['gap'] = df['gap'].apply(lambda x: min(100, x))
    df['gap'] = df['gap'].apply(lambda x: max(0, x))
    # A cause des erreurs d'arrondi, mais toujours positif normalement (TOL=1e-3 peut être un peu trop grand)
    df["closed"] = df["gap"] < 1e-3

    print(df.head())

    # make_graphic(df)
    # return

    ######################
    columns = ["instance"]
    for method in methods:
        if method!='static':
            columns.append(f'objective_{method}')
            columns.append(f'gap_{method}')
        else:
            columns.append('lower_bound_static')
            columns.append('bestlower_bound')
            columns.append('bestUpperBound')

    result_df = pd.DataFrame(columns=columns)

    instance_df_initialied = False
    for instance in df['instance'].unique():
        instance_df = pd.DataFrame({'instance': [instance]})

        for method in methods:
            if method!='static':
                method_df = df[(df['instance'] == instance) & (df['method'] == method)].reset_index(drop=True)
                instance_df[f'objective_{method}'] = method_df['objective'][0]
                instance_df[f'gap_{method}'] = method_df['gap'][0]
            else:
                static_method_df = df[(df['instance'] == instance) & (df['method'] == 'static')].reset_index(drop=True)
                instance_df['lower_bound_static'] = static_method_df['lower_bound'][0]
                instance_df['bestlower_bound'] = df[df['instance'] == instance]['bestlower_bound'].values[0]
                instance_df['bestObjective'] = df[df['instance'] == instance]['bestObjective'].values[0]

        if not instance_df_initialied:
            result_df = instance_df
            instance_df_initialied = True
        else:
            result_df = pd.concat([result_df, instance_df], ignore_index=True)

    # Calcul du poids de la robustesse (inf et sup) différents si l'instance n'a pas été closed
    result_df['PR_inf'] = 100*(result_df['bestlower_bound'] - result_df['lower_bound_static']) / result_df['lower_bound_static']
    result_df['PR_inf'] = result_df['PR_inf'].apply(lambda x: min(100, x))

    result_df['PR_sup'] = 100*(result_df['bestObjective'] - result_df['lower_bound_static']) / result_df['lower_bound_static']
    result_df['PR_sup'] = result_df['PR_sup'].apply(lambda x: min(100, x))
    
    # Define a custom sorting function to sort instances by their numerical value
    def sort_by_instance(instance):
        try:
            return int(instance.split('_')[0].split('/')[-1])
        except ValueError:
            return float('inf')

    # Réorganisation des lignes
    result_df['instance_numeric'] = result_df['instance'].apply(sort_by_instance)
    result_df.sort_values(by='instance_numeric', inplace=True)

    # Réorganisation des colonnes
    final_colums = ['instance', 'PR_inf', 'PR_sup']
    for method in methods:
        if method!='static':
            final_colums.append(f'objective_{method}')
            final_colums.append(f'gap_{method}')
    result_df = result_df[final_colums].reset_index(drop=True)

    # Round the dataframe
    columns_to_round = ['PR_inf', 'PR_sup', 'objective_dualized', 'gap_dualized', 
                    'objective_heuristic', 'gap_heuristic', 'objective_plans_coupants', 
                    'gap_plans_coupants']
    result_df[columns_to_round] = result_df[columns_to_round].round(2).abs()
    # le abs sert juste à enlever les -0.0

    result_df.to_csv('results_.csv', index=False)


def main():
    return


if __name__ == '__main__':
    # main()
    make_results_tab()
