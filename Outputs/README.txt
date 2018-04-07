####################################################################
Readme file describing contents of the 'Outputs' folder
####################################################################
This folder contains output files used in other code files in the 'Code' and 'Validation' folders.

> CSV files
	> Imputed_data_*: contains imputed data- 10 instances.
	> mosq scores: risk scores in the mosquito model.
	> mosq scores: 10-fold cv scores in the mosquito model.
	> FlaviData: covariate data for all flaviviruses
	> mosquito: covariate data for all mosquito species
	> risk scores_mid_NoHuPopDen: Zika risk scores in the primate model.
	> risk scores_10fold_mid_NoHuPopDen: 10-fold Zika risk scores in the primate model.

> RDA files
	> all_models_mosq: all trained mosquito models.
	> all_models_tuned_mid: primate models trained on all 10 imputed datasets.
	> flavi_imputed_data: all imputed datasets.
	> mice_model_best: imputed model distribution parameters. Can be used to generate imputed datasets.
	> pred_matrix_mid_NoHuPopDen: matrix containing primate scores from all 10 primate models.
	> pred_matrix_10fold_mid_NoHuPopDen: matrix containing 10-fold primate scores from all 10 primate models.