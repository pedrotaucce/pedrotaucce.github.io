---
layout: post
title: Building a calibrated gene tree with BEAST 2.7
date: 2024-12-02 11:00:00
description: A step-by-step beginner's guide 
tags: BEAST Dating UltrametricTree MolecularClock
categories: Phylogenetics
thumbnail: assets/img/beast.png
---

Since our blog is called **The BEAST**, it is only fair that we debut our posts by discussing **Bayesian Evolutionary Analysis Sampling Trees**, or simply <a href="https://www.beast2.org">BEAST</a> (<a href = "https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006650">Bouckaert et al. 2019)</a>, one of the most widespread molecular phylogenetic software packages. Although there are many great <a href="https://www.beast2.org/tutorials/">tutorials</a> on the BEAST2 website and even an entire site dedicated to a yearly course taught by some of the developers (see <a href="https://taming-the-beast.org/tutorials/">Taming the BEAST</a>), I believe that a step-by-step guide aimed at new users is always welcome.

In Bayesian phylogenetic analyses, hypotheses (*i.e.*, the tree and its parameters) are generated based on a data matrix, an assumed model of nucleotide evolution, and various priors. We will revisit all of this shortly, but for now, we will need a [molecular matrix](/assets/files/p_rupicola.fasta). For this, we will use the 16S mtDNA matrix from the description of *Pristimantis rupicola* (Anura: Brachycephaloidea), from <a href="https://bioone.org/journals/journal-of-herpetology/volume-54/issue-2/19-114/A-New-Rupicolous-Species-of-the-Pristimantis-conspicillatus-Group-Anura/10.1670/19-114.short">Taucce et al. (2020)</a>. I will assume you are familiar with retrieving sequences from <a href="https://www.ncbi.nlm.nih.gov/genbank/">GenBank</a> and aligning them, but we will cover these topics in future posts.

**BEAST2** interprets XML files, and to build one, we will use the **BEAUti** app, which is part of BEAST software package. First, we will load our molecular matrix by going to <code>File > Import Alignment</code>. BEAUti will prompt us to specify the type of data we have, and we will select "nucleotide" (<a href="#fig1">Fig. 1</a>).

<div class="image-container mt-3">
    <figure id="fig1">
        {% include figure.liquid loading="eager" path="assets/img/posts/2024-10-06/fig1.png" class="img-fluid rounded z-depth-1" zoomable=true %}
         <figcaption style="font-size: 0.9em; text-align: center;"><strong>Figure 1:</strong> BEAUti screen after you have loaded the molecular matrix.</figcaption>
    </figure>
</div>

We will leave the "Tip Dates" panel as it is and move on to the "Site Model" panel (<a href="#fig2">Fig. 2</a>). At this stage, we need to select a nucleotide substitution model for our data. While there are several software packages available for this task, BEAST conveniently includes an inbuilt package that will infer the best-fitting nucleotide substitution model simultaneously with the phylogenetic tree inference: **bModelTest** (<a href="https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-017-0890-6">Bouckaert & Drummond, 2017</a>). If the <code>BEAST Model Test</code> option is not available, you may need to install it by navigating to <code>File > Manage Packages</code>. With that in place, we are now ready to proceed to the next panel, "Clock Model."

<div class="image-container mt-3" style="max-width: 50%; margin: 0 auto;">
    <figure id="fig2">
        {% include figure.liquid loading="eager" path="assets/img/posts/2024-10-06/fig2.png" class="img-fluid rounded z-depth-1" zoomable=true %}
         <figcaption style="font-size: 0.9em; text-align: center;"><strong>Figure 2:</strong> Site Model tab showing the use of the bModelTest package.</figcaption>
    </figure>
</div>

When we open the panel, the strict clock model is already selected. This model assumes that the clock rate evolves equally across the entire tree. However, since our matrix involves interspecific relationships, this assumption is likely not suitable. We could test how "clock-like" our data are, but that will be the subject of a future post. For now, it is advisable to choose the "optimised relaxed clock model".

There are two main ways to date a phylogenetic tree: either by assigning ages to the nodes, usually based on geological events or fossil records, or by specifying a clock rate—*i.e.*, informing BEAST of the nucleotide substitution rate in our matrix. While it is common practice to use a single predetermined clock rate from the literature, we will take a different, more elegant approach. Since it is unlikely that the same clock rate applies to different datasets, estimating the rate is more appropriate. However, to ensure plausible results and aid the analysis in converging, it is recommended to introduce a narrow prior. But how do we do that?

According to <a href="https://doi.org/10.1017/CBO9781139095112">Drummond & Bouckaert (2015)</a>, if "a number of independent rate estimates from other papers on relevant taxa and gene regions" are available, one can fit a log-normal distribution to these values and use them as a prior. Since several estimates for the 16S gene in frogs are available in the literature (see <a href="https://doi.org/10.1111/j.1558-5646.2007.00181.x">Lemmon et al., 2007</a>, <a href="https://doi.org/10.1093/sysbio/syr130">Fouquet et al., 2012</a>, <a href="https://doi.org/10.1016/j.ympev.2017.04.007">Gehara et al., 2017</a>, and <a href="https://doi.org/10.1080/14772000.2022.2102686">Dufresnes et al., 2022</a>), we will follow this approach. All we need is <a href = "https://www.R-project.org/">R</a> (R Core Team, 2024) and a small piece of code:

```r
library(fitdistrplus)
library(stats)

rate <- c(0.006, 0.00555, 0.00277, 0.005)  # Define a vector of rates
ratedist <- fitdist(rate, "lnorm")  # Fit a log-normal distribution to the rate data

mean <- mean(rate)  # Calculate the mean of the rates, it will be useful as a starting point
log.mean <- unname(ratedist$estimate["meanlog"])  # Get the estimated mean value from the fitted lognormal distribution
log.sd <- unname(ratedist$estimate["sdlog"])  # Get the estimated sd value from the fitted lognormal distribution
```

In our case, the mean of the log-normal distribution is -5.374295 and the sd is 0.304072. We will use hese values in the "Priors" panel shortly. For now, let's use the mean of our four rate values, 0.00483, as a starting point to help the Markov chains converge (<a href="#fig3">Fig. 3</a>). **Do not forget to check the "estimate" box**. To do this, navigate to <code>Mode</code> and uncheck the <code>Automatic set clock rate</code> option.

<div class="image-container mt-3" style="max-width: 50%; margin: 0 auto;">
    <figure id="fig3">
        {% include figure.liquid loading="eager" path="assets/img/posts/2024-10-06/fig3.png" class="img-fluid rounded z-depth-1" zoomable=true %}
         <figcaption style="font-size: 0.9em; text-align: center;"><strong>Figure 3:</strong> Clock Model panel showing the chosen model, Optimised Relaxed Clock, and the starting value of the <code>clock.rate</code> prior.</figcaption>
    </figure>
</div>

Moving on to the "Priors" panel, we will see several parameters. We will leave most of them under default, but we will change the tree and the molecular clock priors. The current tree prior, "Yule", has only one parameter to be estimated, the birth rate, assuming no extinction events. A more appropriate prior would be the Birth Death Model (<a href="#fig4">Fig. 4</a>), which does account for extinction.

<div class="image-container mt-3">
    <figure id="fig4">
        {% include figure.liquid loading="eager" path="assets/img/posts/2024-10-06/fig4.png" class="img-fluid rounded z-depth-1" zoomable=true %}
         <figcaption style="font-size: 0.9em; text-align: center;"><strong>Figure 4:</strong> Priors panel, showing the changes in tree and ORCucldMean models.</figcaption>
    </figure>
</div>

Finally, we are going now to use the mean and standard deviation we calculated from our log-normal model. Open the <code>ORCucldMean.c</code> prior and add the meanlog value (in our case, -5.374295) to the M parameter and sdlog (0.304072) to S (<a href="#fig4">Fig. 4</a>). To make sure, you can plot your log-normal model and compare it to the one in BEAUti (<a href="#fig5">Fig. 5</a>):


```r
x <- seq(0, 0.01, length = 100)  # Create a sequence of values from 0 to 0.01 with 100 points
y <- dlnorm(x, meanlog = log.mean, sdlog = log.sd)  # Calculate the density of the log-normal distribution for the x values

# Set up the graphics device to save the plot as a PNG file
png("your_path.png", width = 90, height = 90, units = "mm",
    res = 300, pointsize = 8)

# Plot the density curve of the log-normal distribution
curve(dlnorm(x, meanlog = log.mean, sdlog = log.sd), ylab = "Probability density",
      from = 0.001, to = 0.01, col = "red", xlab = "Clock rate", lwd = 1)

# Add a filled polygon under the density curve
polygon(c(x, rev(x)), c(y, rep(0, length(y))), col = rgb(1, 0, 0, 0.5), border = NA)

dev.off()  # Close the graphics device
```
<div class="image-container mt-3" style="max-width: 50%; margin: 0 auto;">
    <figure id="fig5">
        {% include figure.liquid loading="eager" path="assets/img/posts/2024-10-06/fig5.png" class="img-fluid rounded z-depth-1" zoomable=true %}
        <figcaption style="font-size: 0.9em; text-align: center;"><strong>Figure 5:</strong> Plot of the fitted log-normal model</figcaption>
    </figure>
</div>

We will leave the MCMC panel unchanged. If your dataset is larger, you may want to consider increasing the number of generations, but for our matrix, 10 million generations with sampling every 1,000 generations should be sufficient.

<h3>Results</h3>

We will use <a href="https://github.com/beast-dev/tracer/releases/latest">Tracer</a> (<a href="https://academic.oup.com/sysbio/article/67/5/901/4989127">Rambaut et al., 2018</a>) to assess the convergence of the MCMC output. Ideally, all ESS values should exceed 200. When conducting the final analysis, you should run at least one additional replicate, starting from a different seed, to ensure that the MCMC has fully converged.

To check our analysis, click the "+" button in the top left corner and load the <code>.log</code> file generated by BEAST. As shown in <a href="#fig6">Fig. 6</a>, a chain length of 10,000,000 generations was sufficient for the MCMC output to converge.

<div class="image-container mt-3">
    <figure id="fig6">
        {% include figure.liquid loading="eager" path="assets/img/posts/2024-10-06/fig6.png" class="img-fluid rounded z-depth-1" zoomable=true %}
         <figcaption style="font-size: 0.9em; text-align: center;"><strong>Figure 6:</strong> Tracer output after loading the log file.</figcaption>
    </figure>
</div>

The <code>.trees</code> files generated by BEAST contain multiple trees, representing the posterior distribution from the Bayesian analysis. We will now summarize these into a Maximum Credibility Tree using another BEAST tool, <a href="https://beast2.blogs.auckland.ac.nz/treeannotator/">TreeAnnotator</a>. The only change I will make is regarding node heights, which I will set to "mean heights"—but feel free to use the option that best suits your analysis (<a href="#fig7">Fig. 7</a>).

<div class="image-container mt-3" style="max-width: 50%; margin: 0 auto;">
    <figure id="fig7">
        {% include figure.liquid loading="eager" path="assets/img/posts/2024-10-06/fig7.png" class="img-fluid rounded z-depth-1" zoomable=true %}
         <figcaption style="font-size: 0.9em; text-align: center;"><strong>Figure 7:</strong> TreeAnnotator window with node heights set to "mean heights".</figcaption>
    </figure>
</div>

After adjusting the input to your <code>.trees</code> file and naming the output as you prefer, you can view the results using your favorite tree viewer. I hope your output looks as cool as mine (<a href="#fig8">Fig. 8</a>)!

<div class="image-container mt-3">
    <figure id="fig8">
        {% include figure.liquid loading="eager" path="assets/img/posts/2024-10-06/fig8.png" class="img-fluid rounded z-depth-1" zoomable=true %}
         <figcaption style="font-size: 0.9em; text-align: center;"><strong>Figure 8:</strong> Chronogram of 16S showing the relationships within *Pristimantis* species.</figcaption>
    </figure>
</div>

<h3>References</h3>

Bouckaert, R., & Drummond, A. J. (2017). bModelTest: Bayesian phylogenetic site model averaging and model comparison. *BMC Evolutionary Biology*, **17**(1), 42.

Bouckaert, R., Vaughan, T. G., Barido-Sottani, J., Duchêne, S., Fourment, M., Gavryushkina, A., et al. (2019). BEAST 2.5: An advanced software platform for Bayesian evolutionary analysis. *PLoS Computational Biology*, **15**(4), e1006650.

Dufresnes, C., Mahony, S., Prasad, V. K., Kamei, R. G., Masroor, R., Khan, M. A., … Litvinchuk, S. N. (2022). Shedding light on taxonomic chaos: Diversity and distribution of South Asian skipper frogs (Anura, Dicroglossidae, Euphlyctis). *Systematics and Biodiversity*, **20**(1), 1–25.

Fouquet, A., Noonan, B. P., Rodrigues, M. T., Pech, N., Gilles, A., & Gemmell, N. J. (2012). Multiple Quaternary refugia in the eastern Guiana Shield revealed by comparative phylogeography of 12 frog species. *Systematic Biology*, **61**(3), 461.

Gehara, M., Barth, A., de Oliveira, E. F., Costa, M. A., Haddad, C. F. B., & Vences, M. (2017). Model-based analyses reveal insular population diversification and cryptic frog species in the *Ischnocnema parva* complex in the Atlantic forest of Brazil. *Molecular Phylogenetics and Evolution*, **112**, 68–78.

Lemmon, E. M., Lemmon, A. R., & Cannatella, D. C. (2007). Geological and climatic forces driving speciation in the continentally distributed trilling chorus frogs (*Pseudacris*). *Evolution*, **61**, 2086–2103.

R Core Team (2024). R: A language and environment for statistical computing. *R Foundation for Statistical Computing*, Vienna, Austria.

Taucce, P. P. G., Nascimento, J. S., Trevisan, C. C., Leite, F. S. F., Santana, D. J., Haddad, C. F. B., & Napoli, M. F. (2020). A new rupicolous species of the *Pristimantis conspicillatus* group (Anura: Brachycephaloidea: Craugastoridae) from Central Bahia, Brazil. *Journal of Herpetology*, **54**(2), 245–257.
