cbwsdid_influence_components <- function(model, event_times = NULL){
  if(!inherits(model, "fixest")){
    stop("`model` must be a fixest object returned by `cbwsdid()`.")
  }

  meta <- attr(model, "cbwsdid")
  if(is.null(meta)){
    stop("`model` does not contain `cbwsdid` metadata.")
  }

  components <- meta$if_components
  if(is.null(components)){
    stop(
      "This model does not contain analytic-variance components. ",
      "Refit it with the current version of `cbwsdid()`."
    )
  }

  required.cols <- c("id", "subexp.id", "et", "treated_sa", "b", "R")
  missing.cols <- setdiff(required.cols, names(components))
  if(length(missing.cols) > 0L){
    stop(
      "`if_components` is missing required column(s): ",
      paste(missing.cols, collapse = ", "),
      "."
    )
  }

  available.event.times <- sort(unique(components$et[components$et != -1]))
  if(is.null(event_times)){
    event_times <- available.event.times
  } else {
    event_times <- sort(unique(as.numeric(event_times)))
  }

  missing.event.times <- setdiff(event_times, available.event.times)
  if(length(missing.event.times) > 0L){
    stop(
      "Some requested event times are not available in `if_components`: ",
      paste(missing.event.times, collapse = ", "),
      "."
    )
  }

  list(
    components = components,
    event_times = event_times
  )
}


cbwsdid_influence_cohort <- function(model, event_times = NULL){
  input <- cbwsdid_influence_components(model, event_times = event_times)
  components <- input$components
  event_times <- input$event_times
  cells <- list()

  for(event.time in event_times){
    event.data <- components[
      components$et == event.time &
        !is.na(components$id) &
        is.finite(components$b) &
        is.finite(components$R),
      ,
      drop = FALSE
    ]

    stack.rows <- split(
      seq_len(nrow(event.data)),
      as.character(event.data$subexp.id)
    )

    for(rows in stack.rows){
      stack.data <- event.data[rows, , drop = FALSE]
      treated <- stack.data$treated_sa == 1
      control <- stack.data$treated_sa == 0
      n.treated <- sum(treated)
      control.mass <- sum(stack.data$b[control])

      if(n.treated <= 0L || !is.finite(control.mass) || control.mass <= 0){
        next
      }

      treated.mean <- mean(stack.data$R[treated])
      control.mean <- sum(
        stack.data$b[control] * stack.data$R[control]
      ) / control.mass
      contrast <- treated.mean - control.mean
      contribution <- numeric(nrow(stack.data))
      contribution[treated] <- (
        stack.data$R[treated] - treated.mean
      ) / n.treated
      contribution[control] <- -stack.data$b[control] * (
        stack.data$R[control] - control.mean
      ) / control.mass

      cells[[length(cells) + 1L]] <- tibble::tibble(
        id = stack.data$id,
        subexp.id = stack.data$subexp.id,
        et = event.time,
        treated_sa = stack.data$treated_sa,
        contribution = contribution,
        contrast = contrast,
        n_treated = n.treated
      )
    }
  }

  out <- dplyr::bind_rows(cells)
  if(nrow(out) == 0L){
    stop("No valid treated/control stack cells are available.")
  }

  out
}


cbwsdid_influence <- function(model, event_times = NULL){
  input <- cbwsdid_influence_components(model, event_times = event_times)
  components <- input$components
  event_times <- input$event_times
  cohort.influence <- cbwsdid_influence_cohort(
    model,
    event_times = event_times
  )
  ids <- unique(components$id)
  influence <- matrix(
    0,
    nrow = length(ids),
    ncol = length(event_times),
    dimnames = list(as.character(ids), as.character(event_times))
  )
  theta <- stats::setNames(rep(NA_real_, length(event_times)), event_times)

  for(event.index in seq_along(event_times)){
    event.time <- event_times[event.index]
    event.influence <- cohort.influence[
      cohort.influence$et == event.time,
      ,
      drop = FALSE
    ]
    cell.info <- event.influence %>%
      dplyr::distinct(
        .data$subexp.id,
        .data$n_treated,
        .data$contrast
      )
    n.treated.total <- sum(cell.info$n_treated)
    theta.event <- sum(
      cell.info$n_treated * cell.info$contrast
    ) / n.treated.total
    theta[event.index] <- theta.event

    event.influence$aggregate_contribution <-
      event.influence$n_treated / n.treated.total *
      event.influence$contribution
    treated.rows <- event.influence$treated_sa == 1
    event.influence$aggregate_contribution[treated.rows] <-
      event.influence$aggregate_contribution[treated.rows] +
      (
        event.influence$contrast[treated.rows] - theta.event
      ) / n.treated.total

    unit.contribution <- event.influence %>%
      dplyr::group_by(.data$id) %>%
      dplyr::summarise(
        contribution = sum(.data$aggregate_contribution),
        .groups = "drop"
      )
    unit.position <- match(unit.contribution$id, ids)
    influence[unit.position, event.index] <- unit.contribution$contribution
  }

  attr(influence, "event_times") <- event_times
  attr(influence, "theta") <- theta
  influence
}


cbwsdid_vcov_analytic <- function(model, event_times = NULL){
  influence <- cbwsdid_influence(model, event_times = event_times)
  variance <- crossprod(influence)
  variance <- (variance + t(variance)) / 2

  list(
    vcov = variance,
    std.error = sqrt(pmax(diag(variance), 0)),
    influence = influence,
    event_times = attr(influence, "event_times"),
    estimate = attr(influence, "theta")
  )
}
