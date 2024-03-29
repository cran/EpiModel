#' @title Function to Trigger the End Horizon
#'
#' @inheritParams recovery.net
#' @inherit recovery.net return
#'
#' @details
#' This function triggers the end horizon if a control `end.horizon` exists and
#' its `at` value is equal to the current timestep. The end horizon consists on
#' the removal of a set of modules from the module list.
#'
#' @keywords internal
trigger_end_horizon <- function(dat) {
  end_horizon <- get_control(dat, "end.horizon", override.null.error = TRUE)
  if (is.null(end_horizon) || end_horizon$at != get_current_timestep(dat)) {
    return(dat)
  }

  dat <- remove_modules(dat, end_horizon$modules)
  dat <- set_control(dat, "end.horizon", NULL)

  return(dat)
}

#' Check the format of the end.horizon control
#'
#' @inheritParams recovery.net
#'
#' @return TRUE invisibly on succes, errors otherwise
#'
#' @keywords internal
check_end_horizon_control <- function(dat) {
  end_horizon <- get_control(dat, "end.horizon", override.null.error = TRUE)
  if (is.null(end_horizon))
    invisible(return(TRUE))

  start <- get_control(dat, "start")
  if (!is.numeric(end_horizon$at) || end_horizon$at <= start) {
    stop(
      "The control `end.horizon` `at` field must be greater than the starting",
      " time step: ", start, "\n"
    )
  }

  missing_mods <- setdiff(end_horizon$modules, names(get_modules(dat)))
  if (length(missing_mods) > 0) {
    stop(
      "The control `end.horizon` `modules` field contains modules names not",
      " present in the module list:\n'",
      paste0(missing_mods, collapse = "', '"), "'\n"
    )
  }

   return(invisible(TRUE))
}
