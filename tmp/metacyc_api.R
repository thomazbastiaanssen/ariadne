
library(httr2)
library(xml2)

url <- "https://websvc.biocyc.org/apixml"
fn <- "enzymes-of-reaction"
id <- "META:TRYPSYN-RXN"

req <- request(url) |>
  req_url_query(fn = fn, id = id, detail = "low")

resp <- req_perform(req)
resp_text <- resp_body_xml(resp)

resp_text |>
    xml_find_all(".//Protein") |>
    xml_attr("ID")


fn <- "substrates-of-reaction"

req <- request(url) |>
  req_url_query(fn = fn, id = id, detail = "low")

resp <- req_perform(req)
resp_text <- resp_body_xml(resp)

resp_text |>
    xml_find_all(".//Compound") |>
    xml_attr("ID")


fn <- "enzymes-of-pathway"
id <- "META:GLYCOLYSIS"

req <- request(url) |>
  req_url_query(fn = fn, id = id, detail = "low")

resp <- req_perform(req)
resp_text <- resp_body_xml(resp)

resp_text |>
    xml_find_all(".//Protein") |>
    xml_attr("ID") |>
    na.omit() |>
    as.vector()

fn <- "get-class-direct-subs"

fn <- "pathways-of-enzrxn"
id <- "META:HEX1-RXN"

req <- request(url) |>
  req_url_query(fn = fn, id = id, detail = "low")

resp <- req_perform(req)
resp_text <- resp_body_html(resp)

resp_text |>
    xml_find_all(".//Protein") |>
    xml_attr("ID") |>
    na.omit() |>
    as.vector()

"monomers-of-protein"

fn <- "pathway-components"
pwy <- "META:GLYCOLYSIS"

req <- request(url) |>
  req_url_query(fn = fn, pwy = pwy)

resp <- req_perform(req)
resp_text <- resp_body_xml(resp)

resp_text |>
    xml_find_all(".//Protein") |>
    xml_attr("ID") |>
    na.omit() |>
    as.vector()
