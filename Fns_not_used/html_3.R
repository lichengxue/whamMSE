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
      padding-top: 16px; /* Adjust this value to align the dropdown button */
      padding-bottom: 16px; /* Adjust this value to align the dropdown button */
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
    footer {
    display: flex;
    justify-content: center;
    align-items: center;
    }
    .footer.logo-container {
    display: flex;
    justify-content: center;
    align-items: center;
    }
    .footer img {
    margin-left: 120px;
    }
  </style>
</head>
<body>
  <div class="navbar">
    <a href="index.html"><img src="C:/Users/chengxue.li/whamMSE/icons/home_icon.png" alt="Home" style="vertical-align:middle; margin-right:8px; width:24px; height:24px;" >
    Home</a>
    <div class="dropdown">
      <button class="dropbtn">
      <img src="C:/Users/chengxue.li/whamMSE/icons/vignettes_icon.png" alt="Home" style="vertical-align:middle; margin-right:8px; width:24px; height:24px;">
      Vignettes</button>
      <div class="dropdown-content">'
vignette_links <- sapply(vignettes, function(v) {
  paste0('<a href="', v$url, '">', v$name, '</a>')
})
index_content <- paste0(index_content, paste(vignette_links, collapse = "\n"), '
      </div>
    </div>
    <a href="https://github.com/lichengxue/whamMSE/tree/master/R">
      <img src="C:/Users/chengxue.li/whamMSE/icons/functions_icon.png" alt="Functions" style="vertical-align:middle; margin-right:8px; width:24px; height:24px;" >
      Functions</a>
    <a href="https://github.com/lichengxue/whamMSE">
      <img src="C:/Users/chengxue.li/whamMSE/icons/source_icon.png" alt="Source" style="vertical-align:middle; margin-right:8px; width:24px; height:24px;" >
      Source</a>
    <a href="#news">
      <img src="C:/Users/chengxue.li/whamMSE/icons/news_icon.png" alt="News" style="vertical-align:middle; margin-right:8px; width:24px; height:24px;" >
      News</a>
    <a href="mailto:chengxue.li@noaa.gov">
      <img src="C:/Users/chengxue.li/whamMSE/icons/contact_icon.png" alt="Contact" style="vertical-align:middle; margin-right:8px; width:24px; height:24px;" >
      Contact</a>
    <a href="mailto:chengxue.li@noaa.gov?subject=Bug Report">
      <img src="C:/Users/chengxue.li/whamMSE/icons/bug_report_icon.png" alt="Bug Report" style="vertical-align:middle; margin-right:8px; width:24px; height:24px;" >
      Bug Report</a>
  </div>
  <div class="container">
    <h1>Management Strategy Evaluation Package for WHAM</h1>
    <p> <br>
    The whamMSE package is a comprehensive toolbox designed to conduct management strategy evaluation (MSE)
    for supporting spatial management for complex, heterogeneous populations. The package
    is built upon a state-space age-structured stock assessment model - Woods Hole Assessment Model (WHAM).
    WHAM can incorporate time- and/or age-varying process errors (treated as random effects) on
    population processes, and can also incorporate environmental effects on population processes.
    This package provides tools to simulate and evaluate management strategies (combination of
    data collection methods, stock assessment models, and harvest control rules) by allowing users
    to explore various scenarios and their impacts on estimation and management performance for fish
    stocks with different life-history traits across different regions.
    </p>
    
    <p>
    Capabilities of the whamMSE package:
    </p>
    <ul>
      <ol>1. Different life histories, spatial structures, and movement dynamics</ol>
      <ol>2. Time- and/or age-varying population and fishery processes</ol>
      <ol>3. Different management strategies</ol>
      <ol>4. Comprehensive output analysis tools</ol>
    </ul>
    
    <p> <br>
    For more details about WHAM, please visit: <a href="https://timjmiller.github.io/wham/" target="_blank">wham official website </a>
    </p>
    
  </div>
  <div class="footer">
    <img src="C:/Users/chengxue.li/whamMSE/icons/wham_icon.png" alt="logo 1" style="width:500px; height:300px;">
    <img src="C:/Users/chengxue.li/whamMSE/icons/SPASAM_logo.jpeg" alt="logo 2" style="width:400px; height:400;">
    <img src="C:/Users/chengxue.li/whamMSE/icons/noaa_logo2.png" alt="logo 3" style="width:500px; height:220px;">
  </div>
  
</body>
</html>')


  

writeLines(index_content, "docs/index.html")