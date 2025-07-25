\name{generate_subclustering_report}
\alias{generate_subclustering_report}
\title{Generate Subclustering HTML Report}
\usage{
generate_subclustering_report(csv_path, html_report_path = NULL, model_name = NULL)
}
\arguments{
\item{csv_path}{Character string specifying the path to the subclustering CSV results file}

\item{html_report_path}{Character string specifying the output path for the HTML report. If NULL, will use the same name as csv_path but with .html extension}

\item{model_name}{Character string specifying the model name to display in the report. If NULL, will extract from the CSV filename}
}
\value{
None. This function generates an HTML report file.
}
\description{
Generate a beautiful HTML report for subclustering batch results, showing annotation, reasoning, and top markers in an interactive format.
}
\examples{
\dontrun{
# Generate a subclustering report from CSV results
generate_subclustering_report(
  csv_path = "subclustering_results.csv",
  html_report_path = "subclustering_report.html",
  model_name = "Gemini-2.5-Flash"
)

# Auto-generate HTML filename and model name
generate_subclustering_report("subclustering_results.csv")
}
} 