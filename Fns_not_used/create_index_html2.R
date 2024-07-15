# Render all vignettes
library(rmarkdown)
setwd(getwd())
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
  <title>Package Vignettes</title>
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
    .navbar a {
      float: left;
      display: block;
      color: #f2f2f2;
      text-align: center;
      padding: 14px 16px;
      text-decoration: none;
    }
    .navbar a:hover {
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
      color: white;
      padding: 14px 16px;
      background-color: inherit;
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
    <a href="#">Home</a>
    <div class="dropdown">
      <button class="dropbtn">Vignettes</button>
      <div class="dropdown-content">'
vignette_links <- sapply(vignettes, function(v) {
  paste0('<a href="', v$url, '">', v$name, '</a>')
})
index_content <- paste0(index_content, paste(vignette_links, collapse = "\n"), '
      </div>
    </div>
  </div>
  <div class="container">
    <h1>Package Vignettes</h1>
    <p>Welcome to the documentation for the R package. Here you can find various vignettes that demonstrate how to use the package.</p>
  </div>
</body>
</html>')

writeLines(index_content, "docs/index.html")
