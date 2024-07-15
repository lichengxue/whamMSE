# Render all vignettes
library(rmarkdown)
setwd(here::here())
vignette_files <- list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE)
sapply(vignette_files, render, output_dir = "docs")

# Generate index.html
vignettes <- lapply(vignette_files, function(file) {
  vignette_name <- basename(file)
  vignette_title <- sub("_", " ", sub("\\.Rmd$", "", vignette_name))
  vignette_html <- sub("\\.Rmd$", ".html", vignette_name)
  list(name = vignette_title, url = vignette_html)
})

index_content <- '<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8" />
  <meta name="generator" content="pandoc" />
  <meta http-equiv="X-UA-Compatible" content="IE=EDGE" />
  <title>Management Strategy Evaluation Package for WHAM</title>
  <style>
    body {
      font-family: Arial, sans-serif;
      line-height: 1.6;
      margin: 0;
      padding: 0;
      background-color: #f4f4f4;
    }
    .navbar {
      background-color: #333;
      overflow: hidden;
    }
    .navbar a, .dropbtn {
      float: left;
      display: block;
      color: #f2f2f2;
      text-align: center;
      padding: 14px 16px;
      text-decoration: none;
    }
    .navbar a:hover, .dropdown:hover .dropbtn {
      background-color: #ddd;
      color: black;
    }
    .dropdown {
      float: left;
      overflow: hidden;
    }
    .dropdown .dropbtn {
      cursor: pointer;
      font-size: 16px;
      border: none;
      outline: none;
      background-color: inherit;
      padding-top: 17.5px; /* Adjust this value to align the dropdown button */
      padding-bottom: 17.5px; /* Adjust this value to align the dropdown button */
    }
    .dropdown-content {
      display: none;
      position: absolute;
      background-color: #f9f9f9;
      min-width: 160px;
      box-shadow: 0px 8px 16px 0px rgba(0,0,0,0.2);
      z-index: 1;
    }
    .dropdown-content a {
      float: none;
      color: black;
      padding: 12px 16px;
      text-decoration: none;
      display: block;
      text-align: left;
    }
    .dropdown-content a:hover {
      background-color: #ddd;
    }
    .dropdown:hover .dropdown-content {
      display: block;
    }
    .container {
      max-width: 800px;
      margin: auto;
      background: white;
      padding: 20px;
      box-shadow: 0 0 10px rgba(0, 0, 0, 0.1);
    }
    h1 {
      text-align: center;
      color: #333;
    }
    p {
      text-align: justify;
      margin: 20px 0;
    }
    ul {
      list-style: none;
      padding: 0;
    }
    li {
      margin: 10px 0;
      font-size: 1.1em;
    }
    a {
      text-decoration: none;
      color: #0073e6;
    }
    a:hover {
      text-decoration: underline;
    }
  </style>
</head>
<body>
  <div class="navbar">
    <a href="index.html">Home</a>
    <div class="dropdown">
      <button class="dropbtn">Vignettes</button>
      <div class="dropdown-content">'
vignette_links <- sapply(vignettes, function(v) {
  paste0('<a href="', v$url, '">', v$name, '</a>')
})
index_content <- paste0(index_content, paste(vignette_links, collapse = "\n"), '
      </div>
    </div>
    <a href="#functions">Functions</a>
    <a href="https://github.com/your-username/your-repo-name">Source Code</a>
    <a href="#news">News</a>
    <a href="mailto:your-email@example.com">Contact</a>
    <a href="mailto:your-email@example.com?subject=Bug Report">Bug Report</a>
  </div>
  <div class="container">
    <h1>Management Strategy Evaluation Package for WHAM</h1>
    <p> <br>
    The wham-mse is a comprehensive toolbox designed to conduct management strategy evaluation (MSE)
    for supporting spatial management for complex, heterogeneous populations. The package
    is built upon a state-space age-structured stock assessment model - Woods Hole Assessment Model (WHAM).
    WHAM can incorporate time- and/or age-varying process errors (treated as random effects) on
    population processes, and can also incorporate environmental effects on population processes.
    This package provides tools to simulate and evaluate management strategies (combination of
    data collection methods, stock assessment models, and harvest control rules) by allowing users
    to explore various scenarios and their impacts on estimation and management performance for fish
    stocks with different life-history traits across different regions.
    </p>
  </div>
</body>
</html>')

writeLines(index_content, "docs/index.html")