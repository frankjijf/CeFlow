






    lookup = []
    for v in all_new:
        lbl = v
        if not vars2_bin.empty and v in getattr(vars2_bin, "name", []):
        elif not vars2_nom.empty and v in getattr(vars2_nom, "name", []):
        elif not vars2_ord.empty and v in getattr(vars2_ord, "name", []):
        elif (
            not vars2_cont.empty
            and hasattr(vars2_cont, "name")
            and v in vars2_cont["name"].values
        ):
    lookup.append((v, lbl, v, lbl))
    var_lookup_df = pd.DataFrame(
        lookup, columns=["variable", "label", "orig_var", "orig_label"]
    )

    # 10. Correlation and Excel report
    if dep_var in CE2_Recoded.columns and pd.api.types.is_numeric_dtype(
        CE2_Recoded[dep_var]
    ):
        corr = (
            CE2_Recoded.corr()[dep_var].drop(dep_var).abs().sort_values(ascending=False)
        )
        corr_df = corr.reset_index().rename(
            columns={"index": "variable", dep_var: "correlation"}
        )
    else:
        corr_df = pd.DataFrame(columns=["variable", "correlation"])

    writer = pd.ExcelWriter(
        os.path.join(path_output, "CE2_EDA_report.xlsx"), engine="xlsxwriter"
    )
    if keep_B:
        pd.DataFrame(vars2_bin).assign(new_var=pd.Series(keep_B)).to_excel(
            writer, "Binary", index=False
        )
    if keep_N:
        pd.DataFrame(vars2_nom).assign(new_var=pd.Series(keep_N)).to_excel(
            writer, "Nominal", index=False
        )
    if keep_O:
        pd.DataFrame(vars2_ord).assign(new_var=pd.Series(keep_O)).to_excel(
            writer, "Ordinal", index=False
        )
    if keep_C:
        pd.DataFrame(vars2_cont).assign(new_var=pd.Series(keep_C)).to_excel(
            writer, "Continuous", index=False
        )
    writer.close()













# Macro 3: Variable Reduction macro variables
fast_opt = 'N'  # Fast option will turn off all tests except multivariate regression
# Optional Macro Variables: These all have defaults that can be used
samplesize = 50000  # Sample size to use for variable reduction
redu_weight = 'N'  # Use weights in variable reduction? (Y/N)
sources = 3  # Minimum number of sources to be selected
# Univariate regression option
univ_reg = 'Y'  # Use univariate regression to choose variables? (Y/N)
maxpuni = 0.0001  # Maximum p-value correlation for selecting via univariate regression
# Correlation option
correlation = 'Y'  # Use correlation to choose variables? (Y/N)
corrcut = 0.01  # Minimum correlation between independent variable and dependent variable
# Principal components option
principal = 'Y'  # Use principal components to choose variables? (Y/N)
nprin = 10  # Number of principal components desired
minprin = 0.5  # Minimum factor correlation for selecting via principal component
# Cluster option
cluster = 'Y'  # Use cluster analysis to choose variables? (Y/N)
maxc = 20  # Number of clusters desired
maxratio = 0.5  # Maximum R-squared ratio for selecting variables via clustering
# Linear regression option
regression = 'Y'  # Use linear regression to choose variables? (Y/N)
alphareg = 0.05  # Alpha level for forward selection in linear regression
# Logistic regression option - only applicable if binary dependent variable
logistic = 'Y'  # Use logistic regression to choose variables? (Y/N)
alphalog = 0.05  # Alpha level for forward selection in logistic regression
# Information value option
information = 'Y'  # Use information value to choose variables? (Y/N)
decile = 20  # Number of groups to use when calculating information values
infvcut = 0.01  # Minimum information value for selecting via information value
# Maximum correlation between independent variables option
ind_correlation = 'Y'  # Exclude variables with high correlation to others? (Y/N)
maxcorr = 0.7  # Maximum correlation allowed between independent variables
# Maximum correlation to dependent variable option
ind_dv_corr = 'Y'  # Exclude variables with high correlation to dependent variable? (Y/N)
max_dv_corr = 0.7  # Maximum correlation allowed to dependent variable

# Macro 4: Model selection and tuning
# Optional Macro Variables: These all have defaults that can be used
sel_alpha = 0.05  # Alpha level when selecting/removing variables into the model
includelist = []  # List of variables that have to be in the final model
startlist = []  # List of variables that need to be in the first step
excludelist = []  # List of variables to be excluded from modeling
criteria = 'c'  # Metric to use during variable tuning. Default is c for binary, AdjRsq for other
threshold = 0  # Minimum change in evaluation metric to include variable
SQL_join = 'union'  # Type of join to use between file portions during variable tuning
minimp = 0.01  # Minimum relative importance to keep variable
graph_plot = 'Y'  # Include graphing? (Y/N)

# Macro 5: Final Model selection
# Optional Macro Variables: These all have defaults that can be used
fin_alpha = 0.05  # Alpha level when selecting/removing variables into the model
method = 'stepwise'  # Model build method
fin_num_category = 10  # Maximum number of categories for profiling variables
fin_equal_dist = 'N'  # Use equal distance for dividing variables into groups for profiling (Y/N)









# 2. EDA, profiling and Recode
LogFile = path_output + "02_CE_Var_EDA_Recode_Log_File.log"
LstFile = path_output + "02_CE_Var_EDA_Recode_List_File.lst"
# Placeholder for EDA and recode logic
CE2_Recoded = CE1_Resampled.copy()  # Example, modify as needed

# 3. Variable reduction and ranking
LogFile = path_output + "03_Var_Redu_Log_File.log"
LstFile = path_output + "03_Var_Redu_List_File.lst"
# Placeholder for variable reduction logic
varlist_redu = ['placeholder_var1', 'placeholder_var2']  # Example, replace with actual logic
# Implement actual variable reduction logic here

# 4. Model selection and tuning
LogFile = path_output + "04_Model_Selection_Log_File.log"
LstFile = path_output + "04_Model_Selection_List_File.lst"
# Model selection logic (example: logistic regression for binary DV)
model = LogisticRegression()
X = CE2_Recoded[varlist_redu]
y = CE2_Recoded[dep_var]
model.fit(X, y)

# 5. Final Model build and validation on test sample
LogFile = path_output + "05_Final_Model_Fit_Validation_Log_File.log"
LstFile = path_output + "05_Final_Model_Fit_Validation_LST_File.lst"
# Placeholder for final model validation logic
varlist_final = ['placeholder_final_var1', 'placeholder_final_var2']  # Example, replace with actual logic
X_final = CE2_Recoded[varlist_final]
y_final = CE2_Recoded[dep_var]
final_model = LogisticRegression()
final_model.fit(X_final, y_final)

# Example output
output_predictions = final_model.predict(X_final)

