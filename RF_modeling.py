import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
os.chdir('/public/data/fujingxiang/micro_aging/CMD3_0815/')
import optuna
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import KFold
# define the objective function
def objective_kfold(trial):
    # set parameters space
    n_estimators = trial.suggest_int('n_estimators', 100, 1200)
    #max_features = trial.suggest_categorical('max_features', [1, 'sqrt', 'log2'])
    max_depth = trial.suggest_int('max_depth', 2, 32)
    min_samples_split = trial.suggest_int('min_samples_split', 2, 10)
    min_samples_leaf = trial.suggest_int('min_samples_leaf', 1, 10)
    # 5 fold cross validation
    kf = KFold(n_splits=5, shuffle=True, random_state=123)
    # setup model
    model = RandomForestRegressor(
        #max_features=max_features,
        n_estimators=n_estimators,oob_score=True,bootstrap=True,
        max_depth=max_depth,n_jobs=30,
        min_samples_split=min_samples_split,
        min_samples_leaf=min_samples_leaf

    )

    mae_ls = []
    for train_index, val_index in kf.split(X_train, Y_train):
        x_train_fold, x_val_fold = X_train.iloc[train_index], X_train.iloc[val_index]
        y_train_fold, y_val_fold = Y_train.iloc[train_index], Y_train.iloc[val_index]
        
        # training on the training fold, and validating on the validation fold
        model.fit(x_train_fold, y_train_fold)
        y_pred_val = model.predict(x_val_fold)
        mae = mean_absolute_error(y_val_fold, y_pred_val)
        mae_ls.append(mae)
    mean_score = np.mean(mae_ls)
        
    return mean_score
#plot optimization history
def plot_optimization_history(study):
    plt.figure(figsize=(8, 6))
    plt.plot([t.number for t in study.trials], [t.value for t in study.trials], 'o-', color='black', linewidth=2, markersize=8)
    plt.xlabel('Trial Number', fontsize=14)
    plt.ylabel('Objective Value', fontsize=14)
    plt.tick_params(labelsize=11)
    plt.text(0.7, 0.9, f'Min MAE: {study.best_value:.4f}', fontsize=14, transform=plt.gca().transAxes)
    plt.hlines(study.best_value, 0, len(study.trials), color='darkred', linestyles='--', linewidth=1.5)
    plt.text(study.best_trial.number, study.best_value-0.083, str(study.best_trial.number), ha='center', va='top', fontsize=11, color='darkred')
    plt.title('Optimization History', fontsize=16)
    plt.tight_layout()
    plt.savefig('./ML_modeling/Optimization_History_kfold.pdf', dpi=400,bbox_inches='tight')
# plot microbiota age vs chronological age in training and test set
def plot_age(test,train):
    sns.set(style='white')
    plt.figure(figsize=(8,8))
    sns.regplot(x='age', y='oob_prediction', data=train, 
                scatter_kws={'s':5,'alpha':0.5},line_kws={'label': 'Training'})
    from scipy import stats
    slope, intercept, r_value, p_value, std_err = stats.linregress(train['age'], train['oob_prediction'])
    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(test['age'], test['oob_prediction'])
    plt.text(20, 72, 'Training Fit: y = %0.2f x + %0.2f' % (slope, intercept), fontsize=12)
    plt.text(30, 70, 'R = %0.2f' % r_value, fontsize=12)
    plt.text(20, 68, 'Test Fit: y = %0.2f x + %0.2f' % (slope1, intercept1), fontsize=12)
    plt.text(30, 66, 'R = %0.2f' % r_value1, fontsize=12)
    sns.regplot(x='age', y='oob_prediction', data=test, 
                scatter_kws={'s':5,'alpha':0.5},line_kws={'label': 'Test'},color='darkred')
    plt.xlabel('Chronological Age')
    plt.ylabel('Microbiota Age')
    plt.legend(loc='upper right')
    plt.savefig('./ML_modeling/RF_result_kfold.pdf', dpi=400,bbox_inches='tight')
    plt.close()
#chow test
from statsmodels.formula.api import ols
def chow_test(test_set,train_set):
    test_set['group'] = 'test'
    train_set['group'] = 'train'
    result_all = pd.concat([test_set,train_set],axis=0)
    model_check = ols('oob_prediction ~ age + group + age*group', data=result_all).fit()
    test_res = model_check.summary().tables[1]
    test_res_df = pd.DataFrame(test_res) 
    return test_res_df



if __name__ == '__main__':
    data_input = pd.read_csv("./data_input_for_modeling.csv")
    data_input.set_index('sample_id', inplace=True)
    print("data_input has been created................")


    #split data_input
    X_train,X_test,Y_train,Y_test = train_test_split(data_input.drop(['age'],axis=1),data_input['age'],test_size=0.2,random_state=123)
    print("data_input has been splited,20% test data,80% train data,start to tune hyperparameters............")
    # create a study object and optimize the objective function
    study_kfold = optuna.create_study(direction='minimize')
    # start optimization
    study_kfold.optimize(objective_kfold, n_trials=100,show_progress_bar=True,n_jobs=30)
    # best params and best score
    best_params_kfold = study_kfold.best_params
    best_score_kofld = study_kfold.best_value

    tuned_model = RandomForestRegressor(
            n_estimators=best_params_kfold["n_estimators"],oob_score=True,bootstrap=True,
            max_depth=best_params_kfold["max_depth"],n_jobs=30,
            min_samples_split=best_params_kfold["min_samples_split"],
            min_samples_leaf=best_params_kfold["min_samples_leaf"]
        )
    plot_optimization_history(study_kfold)
    tuned_model.fit(X_train,Y_train)
    tuned_data = pd.DataFrame({'age':Y_train,'oob_prediction':tuned_model.oob_prediction_})
    print("tuned model has been trained, start to predict longitual data and cross disease data............")

    y_pred_test=tuned_model.predict(X_test)
    #combine all data
    result_test = pd.DataFrame({'age':Y_test,'oob_prediction':y_pred_test})
    #plot tuned data and test data microbota age plot
    plot_age(test=result_test,train=tuned_data)
    print("RF result plot has been saved................")
    chow_test_df = chow_test(test_set=result_test,train_set=tuned_data)
    chow_test_df.to_csv('./ML_modeling/chow_test_kfold.csv',index=False,header=False)
    print("chow test result has been saved................")

    # save tuned_data prediction results and result_test,merge into healthy_microbiota_age.csv
    healthy_mage_all = pd.concat([tuned_data,result_test],axis=0)
    print("all result has been saved,start to explain the model............")
    #model explained by shap
    import shap
    explainer = shap.TreeExplainer(tuned_model)
    shap_values = explainer.shap_values(X_test)
    X_test.columns = [i.split('s__')[1] if 's__' in i else i for i in X_test.columns]
    #shape summary plot
    shap.summary_plot(shap_values, X_test, show=False)
    plt.tight_layout()
    plt.rcParams['pdf.fonttype'] = 42
    plt.savefig('./ML_modeling/summary_plot_kfold.pdf',format='pdf',dpi=400)
    plt.close()
    print("summary plot has been saved................")
    #save model
    import dill
    dill.dump_session('/public/data/fujingxiang/micro_aging/CMD3_0815/ML_modeling.pkl')
    print("model has been saved,done................")