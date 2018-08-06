Notes on econometric data for "Costs and consequences of wind turbine wake effects arising from uncoordinated wind energy development" used in the regression analysis and figures.

There are four data files, one for the Texas case in the main text, and three for the other cases in the SI. Key variable descriptions are as follows:

Variables common across all datsets:
PlantID, PlantName - unique EIA identifier and name
C_MW_i - Capacity (MW) for wind farm i (see key below for each case)
MWh_i - Monthly gen (MWh) for wind farm i (also given as Gen_i for some cases)
CF_i - Capacity factor for wind farm i
align_raw, align_spd, align_spd3 - index of monthly wind direction alignment with wind farm orientation, unweighted (raw),
	weighted by wind speed (spd), weighted by wind speed cubed (spd3)
perc_raw, perc_spd, perc_spd3 - fraction of hours in month wind blew +/- 30 degrees of wind farm orientation, unweighted 			(raw), weighted by wind speed (spd), weighted by wind speed cubed (spd3)
CF_i_hat - Predicted capacity factor at windfarm i with wake effects
CF_i_hat_no - Predicted capacity factor at windfarm i without wake effects
CFrat, CFrat_no - Ratio of predicted capacity factor to actual, with and without wake effects
CFdiff, CFdiff_no - Difference in predicted capacity factor to actual, with and without wake effects 
SettlementPointPrice - ERCOT West average hourly real-time market price, weighted by hourly ERCOT wind generation (Texas case only)

Figures plot {CFdiff, CFdiff_no} against month of sample for comparisons of predicted net actuals, and {CF_i_hat, CF_i_hat_no, CF_i} (where i is the downwind farm) against month of sample for comparisons of predicted and actuals. 


Case Specific wind farm i:
Texas: 
R - Roscoe
C - Champion
L - Loraine

Iowa:
WW - Whispering Willow
FC - Franklin County
C - Century

Illinois:
BC - Benton County
ST - Settler's Trail
FR - Fowler Ridge

Kansas:
SH1 - Smoky Hills 1
PR - Post Rock
SH2 - Smoky Hills 2


