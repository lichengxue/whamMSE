# Render all vignettes
library(rmarkdown)
main.dir = "C:/Users/chengxue.li/whamMSE/"
setwd(main.dir)
vignette_files <- list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE)
sapply(vignette_files, render, output_dir = "docs")

news_files <- list.files("news", pattern = "\\.Rmd$", full.names = TRUE)
sapply(news_files, render, output_dir = "docs")

bug_report_files <- list.files("bug_report", pattern = "\\.Rmd$", full.names = TRUE)
sapply(bug_report_files, render, output_dir = "docs")

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
    margin-left: 100px;
    }
    
    .floating-box {
    position: fixed;
    right: 20px;
    top: 100px;
    background-color: #fff;
    border: 1px solid #ccc;
    padding: 10px;
    box-shadow: 0px 0px 10px rgba(0, 0, 0, 0.1);
    z-index: 1000;
    width: 200px;
    }
    
    floating-box p {
    margin: 0;
    font-weight: bold;
    }
    
  </style>
</head>
<body>
  <div class="navbar">
    <a href="index.html"><img src="https://raw.githubusercontent.com/lichengxue/whamMSE/master/icons/home_icon.png" alt="Home" style="vertical-align:middle; margin-right:8px; width:24px; height:24px;" >
    Home</a>
    <div class="dropdown">
      <button class="dropbtn">
      <img src="https://raw.githubusercontent.com/lichengxue/whamMSE/master/icons/vignettes_icon.png" alt="Home" style="vertical-align:middle; margin-right:8px; width:24px; height:24px;">
      Vignettes</button>
      <div class="dropdown-content">'
vignette_links <- sapply(vignettes, function(v) {
  paste0('<a href="', v$url, '">', v$name, '</a>')
})
index_content <- paste0(index_content, paste(vignette_links, collapse = "\n"), '
      </div>
    </div>
    <a href="https://github.com/lichengxue/whamMSE/tree/master/R">
      <img src="https://raw.githubusercontent.com/lichengxue/whamMSE/master/icons/functions_icon.png" alt="Functions" style="vertical-align:middle; margin-right:8px; width:24px; height:24px;" >
      Functions</a>
    <a href="https://github.com/lichengxue/whamMSE">
      <img src="https://raw.githubusercontent.com/lichengxue/whamMSE/master/icons/source_icon.png" alt="Source" style="vertical-align:middle; margin-right:8px; width:24px; height:24px;" >
      Source</a>
    <a href="news_updates.html">
      <img src="https://raw.githubusercontent.com/lichengxue/whamMSE/master/icons/news_icon.png" alt="News" style="vertical-align:middle; margin-right:8px; width:24px; height:24px;" >
      News</a>
    <a href="mailto:chengxue.li@noaa.gov">
      <img src="https://raw.githubusercontent.com/lichengxue/whamMSE/master/icons/contact_icon.png" alt="Contact" style="vertical-align:middle; margin-right:8px; width:24px; height:24px;" >
      Contact</a>
    <a href="bug_report.html">
      <img src="https://raw.githubusercontent.com/lichengxue/whamMSE/master/icons/bug_report_icon.png" alt="Bug Report" style="vertical-align:middle; margin-right:8px; width:24px; height:24px;" >
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
    <strong>Capabilities of the whamMSE package</strong>:
    </p>
    <ul>
      <ol>1. Different fish life histories, spatial structures, and movement dynamics</ol>
      <ol>2. Time- and/or age-varying population and fishery processes</ol>
      <ol>3. Different management strategies</ol>
      <ol>4. Comprehensive output analysis tools</ol>
    </ul>
    
    <p> <br>
    For more details about WHAM, please visit: <a href="https://timjmiller.github.io/wham/" target="_blank">wham official website </a>
    </p>
    
    <p> <strong>NOAA Disclaimer</strong> <br>
    This repository is a scientific product and is not official communication of the 
    National Oceanic and Atmospheric Administration, or the United States Department 
    of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the 
    user assumes responsibility for its use. Any claims against the Department of 
    Commerce or Department of Commerce bureaus stemming from the use of this GitHub 
    project will be governed by all applicable Federal law. Any reference to specific 
    commercial products, processes, or services by service mark, trademark, manufacturer, 
    or otherwise, does not constitute or imply their endorsement, recommendation or 
    favoring by the Department of Commerce. The Department of Commerce seal and logo, 
    or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement 
    of any commercial product or activity by DOC or the United States Government.
    </p>
    
  </div>
  <div class="footer">
    <img src="https://raw.githubusercontent.com/lichengxue/whamMSE/master/icons/wham_icon.png" alt="logo 1" style="width:400px; height:220px;">
    <img src="https://raw.githubusercontent.com/lichengxue/whamMSE/master/icons/SPASAM_logo.jpeg" alt="logo 2" style="width:300px; height:300;">
    <img src="https://raw.githubusercontent.com/lichengxue/whamMSE/master/icons/noaa_logo2.png" alt="logo 3" style="width:400px; height:180px;">
  </div>
  
  <div class="floating-box">
    <p>Developed by: Chengxue Li</p>
    <p>chengxue.li@noaa.gov </p>
    <p>Institution: NOAA Northeast Fisheries Science Center</p>
  </div>
</body>
</html>')


writeLines(index_content, "docs/index.html")