# README
Although the simulation by @lissa_online_2024 has shown that BES seems to perform quite well compared to other methods, one disadvantage is that they used an arbitrary decision criterion to evaluate the pBF. As @hoijtink_tutorial_2019 argues, unlike p-values, Bayes Factors can and should be interpreted on a continuous scale, not as a binary decision. To properly evaluate the general performance of BES, it might thus be better to take the whole pBF distribution into account.

This internship focused on extending the simulation by @lissa_online_2024 to evaluate the performance of BES based on the full distribution of pBF values, aiming to better understand its behaviour across the conditions specified in the original simulation. Each simulation was run with 1000 iterations in RSudio, version 4.4.2 [@RStudio_2024]. To better visualise the pBF distributions, all pBF values were $log_{10}$ transformed. Because the original simulation study mainly focused on sensitivity and specificity, by testing $H_i$ against its complement $H_c$, this internship mainly focused on $pBF_{ic}$ distributions.

Another simulation conducted by @van_wonderen_bayesian_2024 showed that the performance of BES is not always optimal in boundary cases, where neither $H_i$ nor $H_c$ is true. Therefore, this internship focused primarily on these boundary scenarios. In the original simulation, the boundary scenario occurred when $H_i: \rho > 0.1$ was tested against $H_c: \rho \leq 0.1$ and the true effect was $\rho = 0.1$. Given that an estimate can never exactly equal any specific value, by definition both $H_i$ and $H_c$ are not true when the true effect size is exactly 0.1. Since neither hypothesis is true, it was expected that all $log_{10} pBF_ic$ distributions were centered around 0 (corresponding to a pBF of 1 - $log_{10}(1) = 0$). This would indicate neutral evidence for either hypothesis.

After exploring the original $pBF_{ic}$ distributions, this internship extended the original simulation as follows:

- **Simulation 1:**  Testing a wider range of informed hypotheses ($H_i: \rho > -0.8$, $H_i: \rho > 0$ and $H_i: \rho > 0.8$).
- **Simulation 2:** Including other effect size types, namely covariances and regression coefficients ($\beta$). 
- **Simulation 3:** Pooling the results from k studies using RMA and passing the mean correlation estimates and SEs to `bain`, then plot the BFs based on pooled estimates.
- **Simulation 4:** Testing an ordered hypothesis: $H_i: \rho_1 < \rho_2 < \rho_3$. 

All extensions of the original simulations can be found in the forked repository, under [ExtentedSimulationL](https://github.com/lodewijn/BayesianEvidenceSynthesis/tree/master/ExtentedSimulationL).

## Where do I start?

You can load this project in RStudio by opening the file called 'bayesynth.Rproj'.

## Project structure

<!--  You can add rows to this table, using "|" to separate columns.         -->
File            | Description                | Usage         
--------------- | -------------------------- | --------------
README.md       | Description of project     | Human editable
bayesynth.Rproj | Project file               | Loads project 
LICENSE         | User permissions           | Read only     
.worcs          | WORCS metadata YAML        | Read only     
prepare_data.R  | Script to process raw data | Human editable
renv.lock       | Reproducible R environment | Read only     

<!--  You can consider adding the following to this file:                    -->
<!--  * A citation reference for your project                                -->
<!--  * Contact information for questions/comments                           -->
<!--  * How people can offer to contribute to the project                    -->
<!--  * A contributor code of conduct, https://www.contributor-covenant.org/ -->

# Reproducibility

This project uses the Workflow for Open Reproducible Code in Science (WORCS) to
ensure transparency and reproducibility. The workflow is designed to meet the
principles of Open Science throughout a research project. 

To learn how WORCS helps researchers meet the TOP-guidelines and FAIR principles,
read the preprint at https://osf.io/zcvbs/

## WORCS: Advice for authors

* To get started with `worcs`, see the [setup vignette](https://cjvanlissa.github.io/worcs/articles/setup.html)
* For detailed information about the steps of the WORCS workflow, see the [workflow vignette](https://cjvanlissa.github.io/worcs/articles/workflow.html)

## WORCS: Advice for readers

Please refer to the vignette on [reproducing a WORCS project](https://cjvanlissa.github.io/worcs/articles/reproduce.html) for step by step advice.
<!-- If your project deviates from the steps outlined in the vignette on     -->
<!-- reproducing a WORCS project, please provide your own advice for         -->
<!-- readers here.                                                           -->
