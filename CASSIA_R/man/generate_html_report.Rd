\name{generate_html_report}
\alias{generate_html_report}
\title{Generate HTML Evaluation Report}
\usage{
generate_html_report(result_df, gold_col, pred_col, score_col = "score", 
                    reasoning_col = "reasoning", metrics = NULL, 
                    html_report_path = "report.html", model_name = NULL)
}
\arguments{
\item{result_df}{Data frame containing the evaluation results}

\item{gold_col}{Character string specifying the column name for ground truth/gold standard values}

\item{pred_col}{Character string specifying the column name for predicted values}

\item{score_col}{Character string specifying the column name for evaluation scores (default: "score")}

\item{reasoning_col}{Character string specifying the column name for reasoning/explanation text (default: "reasoning")}

\item{metrics}{Named list of pre-calculated metrics (optional, will be calculated if NULL)}

\item{html_report_path}{Character string specifying the output path for HTML report (default: "report.html")}

\item{model_name}{Character string specifying the model name to display in the report (optional)}
}
\value{
None. This function generates an HTML report file.
}
\description{
Generate a comprehensive HTML report for evaluation results with visualizations, metrics, and sample analyses.
}
\examples{
\dontrun{
# Generate evaluation report
generate_html_report(
  result_df = evaluation_results,
  gold_col = "True Cell Type",
  pred_col = "Predicted Cell Type",
  score_col = "evaluation_score",
  reasoning_col = "explanation",
  html_report_path = "evaluation_report.html",
  model_name = "Claude-3.5-Sonnet"
)
}
} 