<h1>DRI R Package</h1>

<p><a href="https://github.com/DUDUNIAO-F/DRI-Rpackage"><img src="https://img.shields.io/badge/R-package-blue.svg" alt="R package"></a></p>

<p><strong>DRI (tsRNA Regulatory Index)</strong> is a generalizable computational framework designed to score tsRNA regulatory network activity from bulk transcriptomes.</p>

<h2>Installation</h2>

<pre><code class="language-r"># Install from GitHub
devtools::install_github("DUDUNIAO-F/DRI-Rpackage")</code></pre>

<h2>Quick Start</h2>

<pre><code class="language-r">library(DRI)

# Use built-in toy data
res <- calculate_DRI(
  expr_mat = expr_toy, 
  tds_genes = c("TG", "TSHR", "PAX8", "SLC5A5", "NKX2-1")
)

# View results
print(res$DRI_score)</code></pre>

<h2>Web Interface</h2>

<p>For non-computational users, a Shiny app is available at: <a href="https://thca.shinyapps.io/dri-pipeline/">https://thca.shinyapps.io/dri-pipeline/</a></p>
