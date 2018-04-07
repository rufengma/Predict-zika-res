####################################################################
Readme file describing contents of the 'Code' folder
####################################################################

> R files
	> BasicStuff*** files: Cover basic exploration of the primate covariate data and different imputation methods <not being used in outputs>
	> bmlml: source codes for the bayesian algorithm
	> Impute_mice: performs imputation through the 'mice' method and generates a model object ../Outputs/mice_model.Rda (the version in which all variables converged is stored as mice_model_best.Rda)
	> MultiModel_MosqFlavi_10Fold: mosquito-virus model- gets in-sample scores, 10-fold CV scores and compares them
	> MultiModel_PrimFlavi_10FoldServer_mid: primate-virus models- gets in-sample and 10-fold scores: server version for multicore
	> MultiModel_PrimFlavi_ValServer: leave-one-out version of MultiModel_PrimFlavi_10FoldServer <not being used now>
	> MultiModel_PrimFlavi_Varimp: variable importance of primate-virus model
	> OutputPlots: expolratory plots for the primate-virus model
	> Spatial_plots: color coded maps based on predictions of the primate-virus model. Static version of the shiny app, no longer used now
	> UniModel_PrimZoonosis: univariate model to predict zoonosis reservoir staus for primates

> RDA files:
	> flavi_small: contains only required covariates for all primates
	> flavi_small_trans: contains required covariates for all primates with appropriate variable transformations