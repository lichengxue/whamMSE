install.packages("usethis")
library(usethis)
library(devtools)

main.dir
setwd(file.path(main.dir,"whamMSE"))
usethis::use_vignette("vignette")

devtools::build_vignettes()
devtools::document()
devtools::build()
devtools::install()

library(whamMSE)

rmarkdown::render("vignettes/Vignettes_Ex1.Rmd",output_dir = "docs")
rmarkdown::render("vignettes/Installation.Rmd",output_dir = "docs")
rmarkdown::render("vignettes/Vignettes_Ex2.Rmd",output_dir = "docs")

library(rmarkdown)

# Set the working directory to your package
setwd(getwd())

# Render all vignettes in the vignettes directory
vignette_files <- list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE)
sapply(vignette_files, rmarkdown::render, output_dir = "docs")

index_content <- "<html>
<head>
  <title>Vignettes Index</title>
</head>
<body>
  <h1>Vignettes Index</h1>
  <ul>
"

# Add links to each vignette
for (vignette_file in vignette_files) {
  vignette_name <- basename(vignette_file)
  vignette_name <- sub("\\.Rmd$", ".html", vignette_name)
  index_content <- paste0(index_content, '    <li><a href="', vignette_name, '">', vignette_name, '</a></li>\n')
}

index_content <- paste0(index_content, "  </ul>
</body>
</html>")

# Write the index.html file to the docs directory
writeLines(index_content, "docs/index.html")

# ---------------------------------

# Create enhanced HTML content
index_content <- "<!DOCTYPE html>
<html>
<head>
  <meta charset='UTF-8'>
  <title>Vignettes Index</title>
  <style>
    body {
      font-family: Arial, sans-serif;
      line-height: 1.6;
      margin: 0;
      padding: 0;
      display: flex;
      justify-content: center;
      align-items: center;
      height: 100vh;
      background-color: #f4f4f4;
    }
    .container {
      max-width: 800px;
      width: 90%;
      background: #fff;
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
  <div class='container'>
    <h1>Vignettes Index</h1>
    <ul>
"

# Add links to each vignette
for (vignette_file in vignette_files) {
  vignette_name <- basename(vignette_file)
  vignette_name_html <- sub("\\.Rmd$", ".html", vignette_name)
  vignette_title <- sub("_", " ", sub("\\.Rmd$", "", vignette_name))
  index_content <- paste0(index_content, '      <li><a href="', vignette_name_html, '">', vignette_title, '</a></li>\n')
}

index_content <- paste0(index_content, "    </ul>
  </div>
</body>
</html>")

# Write the enhanced index.html file to the docs directory
writeLines(index_content, "docs/index.html")