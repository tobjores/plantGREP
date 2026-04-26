# set default pggf options
.pggf_defaults <- structure(
  list(
    out_dir = '.',
    preview_plot = TRUE,
    print_code = TRUE,
    scales = 'fixed',
    space = 'fixed',
    points_below = 10,
    samples = 256,
    trim = TRUE,
    subdensity = FALSE
  ),
  class = 'pggf_options'
)

pggf_defaults <- .pggf_defaults

pggf_config <- function(..., reset_all = FALSE) {
  if (reset_all) {
    pggf_defaults <<- .pggf_defaults
    return(invisible())
  }
  
  arg_list <- list(...)
  
  for (i in seq_along(arg_list)) {
    key = names(arg_list[i])
    value = arg_list[[i]]
    
    if (is.null(key) || key == '') { # `key` not specified
      if (is.null(value)) { # `value` is `NULL` --> ignore
        next
      } else {
        cli::cli_abort(c('Invalid option {.field {value}}', i = 'Options must be specified as {.arg key = value} pairs'))
      }
    } else if (! key %in% names(.pggf_defaults)) { # `key` not a potential option
      cli::cli_abort('Unknown option {.field {key}}') 
    } else {
      fallback <- .pggf_defaults[[key]]
      if (length(value) > length(fallback)) {
        cli::cli_abort(c('Too many values for option {.field {key}}', i = 'Must be a vector of length {length(fallback)}'))
      } else if (is.null(value) || all(is.na(value))) { # `value` is `NULL` or `NA` --> restore default to factory setting
        pggf_defaults[[key]] <<- fallback
      } else if (! all(class(value) == class(fallback))) { # invalid `value`
        cli::cli_abort(c('Invalid value {.field {value}} for option {.field {key}}', i = '{.field {key}} must be of class {.cls {class(fallback)}}'))
      } else { # set option
        pggf_defaults[[key]] <<- value
      }
    }
  }
}

print.pggf_options <- function(options) {
  for (i in seq_along(options)) {
    cat(writeLines(paste0(names(options[i]), ':\t', options[[i]])))
  }
}


# set up a pggfplot
pggfplot <- function(data, filename, x = NULL, y = NULL, facet_col = NULL, facet_row = NULL, facet_wrap = NULL,
                     scales = pggf_defaults$scales, space = pggf_defaults$space, out_dir = pggf_defaults$out_dir, 
                     preview_plot = pggf_defaults$preview_plot, preview_theme = ggplot2::theme_bw(),
                     print_code = pggf_defaults$print_code, xmin = NULL, xmax = NULL,
                     ymin = NULL, ymax = NULL, axis_annotation = NULL) {
  
  # check arguments
  if (! is.character(filename) || length(filename) != 1) {
    cli::cli_abort('{.arg filename} must be a {.cls character} vector of length 1')
  }
  if (! scales %in% c('fixed', 'free_x', 'free_y', 'free', 'square', 'square*')) {
    cli::cli_abort('{.arg scales} must be one of "fixed", "free_x", "free_y", "free", "square", or "square*"')
  }
  if (! space %in% c('fixed', 'free_x', 'free_y', 'free')) {
    cli::cli_abort('{.arg space} must be one of "fixed", "free_x", "free_y", or "free"')
  }
  
  # convert `x` and `y` to strings
  x_str <- if (is.symbol(substitute(x))) {deparse(substitute(x))} else if (is.null(x)) {x_str <- NULL} else {as.character(x)}
  y_str <- if (is.symbol(substitute(y))) {deparse(substitute(y))} else if (is.null(y)) {y_str <- NULL} else {as.character(y)}
  
  # check if at least one of `x` or `y` is given
  if (is.null(x_str) && is.null(y_str)) {
    cli::cli_abort('At least one of {.arg x} or {.arg y} must be given')
  }
  
  # prepare data
  data <- data |>
    tibble::as_tibble() |>
    dplyr::ungroup()
  
  # get facet info
  if (
    (missing(facet_col) || .pggf_is_null(substitute(facet_col))) && 
    (missing(facet_row) || .pggf_is_null(substitute(facet_row)))
  ) { # neither `facet_col` nor `facet_row` given
    if (! missing(facet_wrap) && ! .pggf_is_null(substitute(facet_wrap))) { # `facet_wrap` given
      wrap <- data |> dplyr::pull({{facet_wrap}}) |> unique() |> sort()
      facets <- tibble::tibble(wrap = wrap) |> tibble::rowid_to_column('facet_id')
      data <- data |> dplyr::inner_join(facets, dplyr::join_by({{facet_wrap}} == wrap))
    } else { # all three missing
      wrap <- NULL
      facets <- NULL
    }
    cols <- NULL
    rows <- NULL
  } else { # `facet_col` and/or `facet_row` given
    if (! missing(facet_wrap) && ! .pggf_is_null(substitute(facet_wrap))) { # `facet_wrap` also given
      cli::cli_abort('{.arg facet_wrap} cannot be specified together with {.arg facet_col} or {.arg facet_row}')
    } else { # `facet_wrap` not given
      wrap <- NULL
      if (! missing(facet_col) && ! .pggf_is_null(substitute(facet_col))) {
        cols <- data |> dplyr::pull({{facet_col}}) |> unique() |> sort()
      } else {
        cols <- NULL
      }
      if (! missing(facet_row) && ! .pggf_is_null(substitute(facet_row))) {
        rows <- data |> dplyr::pull({{facet_row}}) |> unique() |> sort()
      } else {
        rows <- NULL
      }
      facets <- tidyr::expand_grid('row' = rows, 'col' = cols) |> tibble::rowid_to_column('facet_id')
      if (is.null(rows)) { # only `facet_col` given
        data <- data |> dplyr::inner_join(facets, dplyr::join_by({{facet_col}} == col))
      } else if (is.null(cols)) { # only `facet_row` given
        data <- data |> dplyr::inner_join(facets, dplyr::join_by({{facet_row}} == row))
      } else { # both given
        data <- data |> dplyr::inner_join(facets, dplyr::join_by({{facet_col}} == col, {{facet_row}} == row))
      }
    }
  }
  
  # rename `x` and `y` columns
  data <- data |>
    dplyr::rename('x' = {{x}}, 'y' = {{y}})
  
  # convert non-numerical `x` and `y` data
  if (! 'x' %in% names(data) || is.numeric(data$x)) {
    numeric_x <- TRUE
  } else {
    if (scales %in% c('square', 'square*')) {
      cli::cli_abort('Option {.arg scales = "square(*)"} is not compatible with non-numeric data')
    }
    
    numeric_x <- FALSE
    
    data$x <- as.ordered(data$x)
    
    if (is.null(facets) || ! is.null(wrap) || scales %in% c('fixed', 'free_y') || is.null(cols)) {
      data <- data |>
        dplyr::mutate(
          xticklabels = as.character(x),
          x = as.numeric(droplevels(x))
        )
    } else {
      data <- data |>
        dplyr::group_by({{facet_col}}) |>
        dplyr::mutate(
          xticklabels = as.character(x),
          x = as.numeric(droplevels(x))
        ) |>
        dplyr::ungroup()
    }
  }
  if (! 'y' %in% names(data) || is.numeric(data$y)) {
    numeric_y <- TRUE
  } else {
    if (scales %in% c('square', 'square*')) {
      cli::cli_abort('Option {.arg scales = "square(*)"} is not compatible with non-numeric data')
    }
    
    numeric_y <- FALSE
    
    data$y <- as.ordered(data$y)
    
    if (is.null(facets) || ! is.null(wrap) || scales %in% c('fixed', 'free_x') || is.null(rows)) {
      data <- data |>
        dplyr::mutate(
          yticklabels = as.character(y),
          y = as.numeric(droplevels(y))
        )
    } else {
      data <- data |>
        dplyr::group_by({{facet_row}}) |>
        dplyr::mutate(
          yticklabels = as.character(y),
          y = as.numeric(droplevels(y))
        ) |>
        dplyr::ungroup()
    }
  }
  
  # store supplied min/max values
  if (! is.null(xmin)) {
    if (! numeric_x) {
      cli::cli_abort('{.arg xmin} can only be supplied for numeric x values')
    } else {
      if (length(xmin) %in% c(1, length(cols))) {
        xmin <- tibble::tibble('col' = cols, 'x' = xmin)
      } else {
        cli::cli_abort(
          paste0(
            '{.arg xmin} must be a {.cls numeric} vector of length 1',
            ifelse(is.null(cols), '', ' or {length(cols)} (number of facet columns)')
          )
        )
      }
    }
  }
  if (! is.null(xmax)) {
    if (! numeric_x) {
      cli::cli_abort('{.arg xmax} can only be supplied for numeric x values')
    } else {
      if (length(xmax) %in% c(1, length(cols))) {
        xmax <- tibble::tibble('col' = cols, 'x' = xmax)
      } else {
        cli::cli_abort(
          paste0(
            '{.arg xmax} must be a {.cls numeric} vector of length 1',
            ifelse(is.null(cols), '', ' or {length(cols)} (number of facet columns)')
          )
        )
      }
    }
  }
  extra_x <- dplyr::bind_rows(xmin, xmax)
  if (nrow(extra_x) == 0) {extra_x <- NULL}
  
  if (! is.null(ymin)) {
    if (! numeric_y) {
      cli::cli_abort('{.arg ymin} can only be supplied for numeric y values')
    } else {
      if (length(ymin) %in% c(1, length(rows))) {
        ymin <- tibble::tibble('row' = rows, 'y' = ymin)
      } else {
        cli::cli_abort(
          paste0(
            '{.arg ymin} must be a {.cls numeric} vector of length 1',
            ifelse(is.null(rows), '', ' or {length(rows)} (number of facet rows)')
          )
        )
      }
    }
  }
  if (! is.null(ymax)) {
    if (! numeric_y) {
      cli::cli_abort('{.arg ymax} can only be supplied for numeric y values')
    } else {
      if (length(ymax) %in% c(1, length(rows))) {
        ymax <- tibble::tibble('row' = rows, 'y' = ymax)
      } else {
        cli::cli_abort(
          paste0(
            '{.arg ymax} must be a {.cls numeric} vector of length 1',
            ifelse(is.null(rows), '', ' or {length(rows)} (number of facet rows)')
          )
        )
      }
    }
  }
  extra_y <- dplyr::bind_rows(ymin, ymax)
  if (nrow(extra_y) == 0) {extra_y <- NULL}
  
  if (! is.null(extra_x) || ! is.null(extra_y)) {
    if (is.null(facets)) {
      extra_coords <- dplyr::bind_cols(extra_x, extra_y)
      if (nrow(extra_coords) == 0) {extra_coords <- NULL}
    } else {
      if (is.null(extra_x)) {
        extra_coords <- facets
      } else if (is.null(cols)) {
        extra_coords <- facets |>
          dplyr::cross_join(extra_x)
      } else {
        extra_coords <- facets |>
          dplyr::inner_join(extra_x, by = 'col', relationship = 'many-to-many')
      }
      if (! is.null(extra_y)) {
        if (is.null(rows) || is.null(extra_y)) {
          extra_coords <- extra_coords |>
            dplyr::cross_join(extra_y)
        } else {
          extra_coords <- extra_coords |>
            dplyr::inner_join(extra_y, by = 'row', relationship = 'many-to-many')
        }
      }
      if ('col' %in% names(extra_coords)) {
        extra_coords <- extra_coords |>
          dplyr::rename({{facet_col}} := col)
      }
      if ('row' %in% names(extra_coords)) {
        extra_coords <- extra_coords |>
          dplyr::rename({{facet_row}} := row)
      }
    }
  } else {
    extra_coords <- NULL
  }
  
  # add axis annotation
  if (! is.null(axis_annotation)) {
    if (is.list(axis_annotation) && ! is.null(names(axis_annotation)) && all(names(axis_annotation) != '') && all(sapply(axis_annotation, length) %in% c(1, nrow(facets)))) {
      if (! is.null(facets)) {
        axis_annotation <- facets |>
          dplyr::select(facet_id) |>
          bind_cols(axis_annotation)
      } else {
        axis_annotation <- as_tibble(axis_annotation)
      }
    } else {
      cli::cli_abort(
        paste0(
          '{.arg axis_annotation} must be a named {.cls list} of length 1',
          ifelse(is.null(facets), '', ' or {nrow(facets)} (number of facets)')
        )
      )
    }
  }
  
  # create a preview plot
  if (preview_plot) {
    plot <- list(preview_theme)
    if (! is.null(facets)) {
      if (! is.null(wrap)) {
        plot_facet_wrap <- .pggf_sym_or_str(substitute(facet_wrap))
        plot <- append(plot, ggplot2::facet_wrap(ggplot2::vars({{plot_facet_wrap}})))
      } else {
        plot_facet_col <- .pggf_sym_or_str(substitute(facet_col))
        plot_facet_row <- .pggf_sym_or_str(substitute(facet_row))
        if (scales %in% c('square', 'square*')) {
          plot_scales <- 'fixed'
        } else {
          plot_scales <- scales
        }
        plot <- append(
          plot,
          ggplot2::facet_grid(
            cols = ggplot2::vars({{plot_facet_col}}),
            rows = ggplot2::vars({{plot_facet_row}}),
            scales = plot_scales,
            space = space
          )
        )
      }
    }
  } else {
    plot <- NULL
  }
  
  # generate LaTeX code
  if (print_code) {
    if (! is.null(wrap)) {
      code_wrap <- '\twrap cols = 3, % adjust to change number of columns'
    } else {
      code_wrap <- NULL
    }
    code <- c('\\begin{pggfplot}[%', code_wrap, paste0('\t% <pggfplot options>\n]{', filename, '}\n'), '\\end{pggfplot}')
  } else {
    code <- NULL
  }
  
  # set up `pggf` object
  pggf <- structure(
    list(
      data = data,
      facets = facets,
      facet_col = rlang::enquo(facet_col),
      facet_row = rlang::enquo(facet_row),
      facet_wrap = rlang::enquo(facet_wrap),
      scales = scales,
      space = space,
      error = list('pggf_x_err_minus' = 0, 'pggf_x_err_plus' = 0, 'pggf_y_err_minus' = 0, 'pggf_y_err_plus' = 0),
      extra_coords = extra_coords,
      numeric_x = numeric_x,
      numeric_y = numeric_y,
      x_str = x_str,
      y_str = y_str,
      plot = plot,
      code = code,
      filename = filename,
      out_dir = out_dir,
      files = NULL,
      stats = NULL,
      axes_extra = NULL,
      axis_annotation = axis_annotation
    ),
    class = 'pggf'
  )
  
  # call `pggf` object
  pggf
}


# "print" `pggf` objects
print.pggf <- function(pggf) {
  # make sure `pggf` is a `pggf` class object
  .assert_pggf(pggf, check_plots = FALSE)
  
  # check if there are any plot files
  if (is.null(pggf$files)) {
    cli::cli_abort('No plots added to {.fn pggfplot} yet.')
  } else {
    # export plot files
    for (file in names(pggf$files)) {
      readr::write_file(pggf$files[file], file.path(pggf$out_dir, file))
    }
  }
  
  # export stats
  if (! is.null(pggf$stats)) {
    readr::write_file(pggf$stats, file.path(pggf$out_dir, paste0(pggf$filename, '_stats.tsv')))
  }
  
  # export facet info
  if (! is.null(pggf$facets)) {
    readr::write_tsv(pggf$facets, file.path(pggf$out_dir, paste0(pggf$filename, '_facets.tsv')))
  }
  
  # calculate and export axis limits
  axis_limits <- .pggf_get_axes(pggf)
  
  # adjust axis limits for ridgeline plots
  if (paste0(pggf$filename, '_ridgeline.tsv') %in% names(pggf$files)) {
    ridge_axis <- gsub('(^|.*\\t)density_(.)(\\t|\\n).*', '\\2', pggf$files[paste0(pggf$filename, '_ridgeline.tsv')])
    
    axis_limits <- axis_limits |>
      dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
      dplyr::mutate(
        .pggf.line_id = dplyr::cur_group_id()
      ) |>
      dplyr::ungroup()
    
    if ('facet_id' %in% names(axis_limits)) {
      axis_helpers <- pggf$data |>
        dplyr::mutate(
          facet_id = as.character(facet_id)
        ) |>
        dplyr::inner_join(
          axis_limits |>
            dplyr::select(facet_id, .pggf.line_id) |>
            tidyr::separate_longer_delim(
              facet_id,
              delim = ','
            ),
          by = 'facet_id'
        )
    } else {
      axis_helpers <- pggf$data |>
        dplyr::mutate(
          .pggf.line_id = 1
        )
    }
    
    axis_helpers <- axis_helpers |>
      dplyr::group_by(.pggf.line_id) |>
      dplyr::summarise(
        .pggf.min = min(.data[[ridge_axis]]),
        .pggf.max = max(.data[[ridge_axis]]),
        .pggf.dens = max(density[.data[[ridge_axis]] == max(.data[[ridge_axis]])]),
        .groups = 'drop'
      ) |>
      dplyr::ungroup()
    
    axis_limits <- axis_limits |>
      dplyr::inner_join(axis_helpers, by = '.pggf.line_id') |>
      dplyr::mutate(
        !! paste0(ridge_axis, 'min') := .pggf.min,
        !! paste0(ridge_axis, 'max') := paste0(.pggf.max, '+', .pggf.dens, '*\\pggfridgescale')
      ) |>
      dplyr::select(-tidyselect::starts_with('.pggf.'))
  } else {
    ridge_axis <- NULL
  }
  
  readr::write_tsv(axis_limits, file.path(pggf$out_dir, paste0(pggf$filename, '_axes.tsv')))
  
  # show preview plot
  if (! is.null(pggf$plot)) {
    # set ridgeline scale to 1 for the preview plot
    if (! is.null(ridge_axis)) {
      axis_limits <- axis_limits |>
        dplyr::rowwise() |>
        dplyr::mutate(
          dplyr::across(
            tidyselect::all_of(paste0(ridge_axis, c('min', 'max'))),
            ~ eval(parse(text = gsub('\\\\pggfridgescale', '1', .x)))
          )
        ) |>
        dplyr::ungroup()
    }
    
    # update preview plot axis limits and labels
    if (pggf$scales %in% c('free_x', 'free') && nrow(axis_limits > 1) && 'col' %in% names(pggf$facets)) {
      axis_limits <- axis_limits |>
        tidyr::separate_longer_delim(facet_id, ',') |>
        dplyr::mutate(facet_id = as.numeric(facet_id)) |>
        dplyr::inner_join(pggf$facets, by = 'facet_id') |>
        dplyr::rename(!! pggf$facet_col := col)
    }
    pggf$plot <- append(
      pggf$plot,
      c(
        ggplot2::expand_limits(axis_limits |> dplyr::select(-ymin, -ymax))
      )
    )
    
    if (! pggf$numeric_x) {
      # replace x with xticklabels
      x_order <- pggf$data |>
        dplyr::group_by(xticklabels) |>
        dplyr::summarise(x = max(x), .groups = 'drop') |>
        dplyr::arrange(x) |>
        dplyr::pull(xticklabels)
      
      pggf$data <- pggf$data |>
        dplyr::mutate(
          x.orig = x,
          x = ordered(xticklabels, levels = x_order)
        ) |>
        dplyr::select(-xticklabels)
      
      pggf$plot <- append(
        pggf$plot,
        c(
          ggplot2::scale_x_discrete(expand = c(0, 0))
        )
      )
    }

    if (pggf$scales %in% c('free_y', 'free') && nrow(axis_limits > 1) && 'row' %in% names(pggf$facets)) {
      if (! 'row' %in% names(axis_limits)) {
        axis_limits <- axis_limits |>
          tidyr::separate_longer_delim(facet_id, ',') |>
          dplyr::mutate(facet_id = as.numeric(facet_id)) |>
          dplyr::inner_join(pggf$facets, by = 'facet_id')
      }
      axis_limits <- axis_limits |>
        dplyr::rename(!! pggf$facet_row := row)
    }
    pggf$plot <- append(
      pggf$plot,
      c(
        ggplot2::expand_limits(axis_limits |> dplyr::select(-xmin, -xmax))
      )
    ) 

    if (! pggf$numeric_y) {
      # replace y with yticklabels
      y_order <- pggf$data |>
        dplyr::group_by(yticklabels) |>
        dplyr::summarise(y = max(y), .groups = 'drop') |>
        dplyr::arrange(y) |>
        dplyr::pull(yticklabels)
      
      pggf$data <- pggf$data |>
        dplyr::mutate(
          y.orig = y,
          y = ordered(yticklabels, levels = y_order)
        ) |>
        dplyr::select(-yticklabels)
      
      pggf$plot <- append(
        pggf$plot,
        c(
          ggplot2::scale_y_discrete(expand = c(0, 0))
        )
      )
    }
    
    plot <- ggplot2::ggplot(pggf$data, ggplot2::aes(x = x, y = y)) +
      ggplot2::labs(
        x = pggf$x_str,
        y = pggf$y_str
      ) +
      pggf$plot
    
    print(plot)
  }
  
  # print LaTeX code
  if (! is.null(pggf$code)) {
    # add x and y labels
    pggf$code <- append(pggf$code, paste0(c('\txlabel = ', '\tylabel = '), c(pggf$x_str, pggf$y_str), ','), after = 1)
    cat(writeLines(pggf$code))
  }
  
  # return `pggf` object invisibly
  invisible(pggf)
}


# calculate axis limits
.pggf_get_axes <- function(pggf) {
  # make sure `pggf` is a `pggf` class object
  .assert_pggf(pggf, dims = 2, check_plots = FALSE)
  
  # prepare data
  if (paste0(pggf$filename, '_hexbin.tsv') %in% names(pggf$files)) {
    axis_limits <- pggf$files[[paste0(pggf$filename, '_hexbin.tsv')]] |>
      readr::read_tsv()
  } else {
    axis_limits <- pggf$data
  }
  
  axis_limits <- axis_limits |>
    dplyr::bind_cols(pggf$error) |>
    dplyr::bind_rows(pggf$extra_coords) |>
    dplyr::mutate(
      dplyr::across(dplyr::matches('pggf_._err_.*'), ~ tidyr::replace_na(.x, 0))
    )
  
  # protect ticklabels with commas
  axis_limits <- axis_limits |>
    dplyr::mutate(
      dplyr::across(tidyselect::ends_with('ticklabels'), ~ ifelse(grepl(',', .x, fixed = TRUE), paste0('{', .x, '}'), .x))
    )
  
  # add all possible facets
  if (! is.null(pggf$facets) && ! 'wrap' %in% names(pggf$facets)) {
    if (! 'row' %in% names(pggf$facets)) { # only `facet_col` given
      has_cols <- TRUE
      has_rows <- FALSE
      axis_limits <- axis_limits |> dplyr::full_join(pggf$facets, dplyr::join_by(facet_id, !! pggf$facet_col == col))
    } else if (! 'col' %in% names(pggf$facets)) { # only `facet_row` given
      has_cols <- FALSE
      has_rows <- TRUE
      axis_limits <- axis_limits |> dplyr::full_join(pggf$facets, dplyr::join_by(facet_id, !! pggf$facet_row == row))
    } else { # both given
      has_cols <- TRUE
      has_rows <- TRUE
      axis_limits <- axis_limits |> dplyr::full_join(pggf$facets, dplyr::join_by(facet_id, !! pggf$facet_col == col, !! pggf$facet_row == row))
    }
  } else {
    has_cols <- FALSE
    has_rows <- FALSE
  }
  
  # get x limits (and labels for non-numeric data)
  if (! is.null(pggf$facets) && has_cols && pggf$scales %in% c('free_x', 'free')) {
    axis_limits <- axis_limits |>
      dplyr::group_by(!! pggf$facet_col)
  }
  
  axis_limits <- axis_limits |>
    dplyr::mutate(
      dplyr::across(x, list('min' = ~ min(.x - pggf_x_err_minus, na.rm = TRUE), 'max' = ~ max(.x + pggf_x_err_plus, na.rm = TRUE)), .names = '{.col}{.fn}'),
      dplyr::across(tidyselect::any_of('xticklabels'), ~ paste0(unique(.x[order(x)]), collapse = ','))
    ) |>
    dplyr::ungroup()
  
  # get y limits (and labels for non-numeric data)
  if (! is.null(pggf$facets) && has_rows && pggf$scales %in% c('free_y', 'free')) {
    axis_limits <- axis_limits |>
      dplyr::group_by(!! pggf$facet_row)
  }
  
  axis_limits <- axis_limits |>
    dplyr::mutate(
      dplyr::across(y, list('min' = ~ min(.x - pggf_y_err_minus, na.rm = TRUE), 'max' = ~ max(.x + pggf_y_err_plus, na.rm = TRUE)), .names = '{.col}{.fn}'),
      dplyr::across(tidyselect::any_of('yticklabels'), ~ paste0(unique(.x[order(y)]), collapse = ','))
    ) |>
    dplyr::ungroup()
  
  # select distinct rows
  axis_limits <- axis_limits |>
    dplyr::distinct(xmin, xmax, ymin, ymax, dplyr::pick(tidyselect::any_of(c('xticklabels', 'yticklabels', 'facet_id'))))
  
  # get x and y ticks (only for non-numeric data)
  if (! pggf$numeric_x) {
    axis_limits <- axis_limits |>
      dplyr::rowwise() |>
      dplyr::mutate(xtick = paste0(seq(xmin, xmax), collapse = ',')) |>
      dplyr::ungroup()
  }
  if (! pggf$numeric_y) {
    axis_limits <- axis_limits |>
      dplyr::rowwise() |>
      dplyr::mutate(ytick = paste0(seq(ymin, ymax), collapse = ',')) |>
      dplyr::ungroup()
  }
  
  # add extra columns to axis limits
  if (! is.null(pggf$axes_extra)) {
    if ('facet_id' %in% names(axis_limits)) {
      axis_limits <- axis_limits |>
        dplyr::inner_join(
          pggf$axes_extra,
          by = 'facet_id'
        )
    } else {
      axis_limits <- axis_limits |>
        dplyr::bind_cols(pggf$axes_extra)
    }
  }
  
  # add axis annotation to axis limits
  if (! is.null(pggf$axis_annotation)) {
    if ('facet_id' %in% names(axis_limits)) {
      axis_limits <- axis_limits |>
        dplyr::full_join(
          pggf$axis_annotation,
          by = 'facet_id'
        )
    } else {
      axis_limits <- axis_limits |>
        dplyr::bind_cols(pggf$axis_annotation)
    }
  }
  
  # enlarge axis limits for hexbin plots
  if (paste0(pggf$filename, '_hexbin.tsv') %in% names(pggf$files)) {
    axis_limits <- axis_limits |>
      dplyr::mutate(
        .pggf.extrax = 0.5 * .pggf.hex.width,
        xmin = xmin - .pggf.extrax,
        xmax = xmax + .pggf.extrax,
        .pggf.extray = 2/3 * (ymax - ymin) * (.pggf.dim.y / (.pggf.dim.y - 1) - 1),
        ymin = ymin - .pggf.extray,
        ymax = ymax + .pggf.extray,
        `hexbin width` = .pggf.hex.width,
        `hexbin shape` = .pggf.hex.shape
      ) |>
      dplyr::select(-tidyselect::starts_with('.pggf.'))
  }
  
  # enforce equal x and y axis limits for `scales = "square"`
  if (pggf$scales == 'square') {
    axis_limits <- axis_limits |>
      dplyr::mutate(
        xmin = min(xmin, ymin, na.rm = TRUE),
        xmax = max(xmax, ymax, na.rm = TRUE),
        ymin = xmin,
        ymax = xmax,
      )
  }
  
  # enforce equal x and y axis limits for `scales = "square*"`
  if (pggf$scales == 'square*') {
    axis_limits <- axis_limits |>
      dplyr::mutate(
        xmin = min(xmin, -ymax, na.rm = TRUE),
        xmax = max(xmax, -ymin, na.rm = TRUE),
        ymin = -xmax,
        ymax = -xmin,
      )
  }
  
  # enlarge axis limits for non-numeric coordinates
  if (! pggf$numeric_x) {
    axis_limits <- axis_limits |>
      dplyr::mutate(
        xmin = xmin - 0.5,
        xmax = xmax + 0.5,
      )
  }
  if (! pggf$numeric_y) {
    axis_limits <- axis_limits |>
      dplyr::mutate(
        ymin = ymin - 0.5,
        ymax = ymax + 0.5,
      )
  }
  
  axis_limits <- axis_limits |>
    dplyr::select(-tidyselect::ends_with('diff'))
  
  # scale axes relative to their limits
  if (has_cols && pggf$space %in% c('free_x', 'free') && pggf$scales %in% c('free_x', 'free')) {
    axis_limits <- axis_limits |>
      dplyr::mutate(
        `scale facet width` = (xmax - xmin) / sum(xmax - xmin) * dplyr::n()
      )
  }
  if (has_rows && pggf$space %in% c('free_y', 'free') && pggf$scales %in% c('free_y', 'free')) {
    axis_limits <- axis_limits |>
      dplyr::mutate(
        `scale facet height` = (ymax - ymin) / sum(ymax - ymin) * dplyr::n()
      )
  }
  
  # collapse identical rows
  axis_limits <- axis_limits |>
    dplyr::group_by(dplyr::pick(-tidyselect::any_of('facet_id'))) |>
    dplyr::summarise(
      dplyr::across(tidyselect::any_of('facet_id'), ~ paste0(sort(unique(.x)), collapse = ',')),
      .groups = 'drop'
    ) |>
    dplyr::arrange(dplyr::pick(tidyselect::any_of('facet_id')))
  
  # return axis limits
  return(axis_limits)
}


# test if an object is a valid `pggf` class object
.assert_pggf <- function(pggf, dims = 0, check_plots = TRUE) {
  if (! inherits(pggf, 'pggf')) {
    cli::cli_abort('{.arg pggf} must be a {.cls pggf} object')
  }
  
  if (dims == 2 && (is.null(pggf$x_str) || is.null(pggf$y_str))) {
    fn <- as.character(sys.call(sys.parent(1)))[1]
    cli::cli_abort('Both {.arg x} and {.arg y} must be given for {.fn {fn}}')
  } else if (dims == 1 && ! is.null(pggf$x_str) && ! is.null(pggf$y_str)) {
    fn <- as.character(sys.call(sys.parent(1)))[1]
    cli::cli_abort('Only {.arg x} or {.arg y} can be given for {.fn {fn}}')
  }
  
  # test if `pggf_density`, `pggf_ridgeline` or `pggf_histogram` has been called before
  if (check_plots && ! is.null(pggf$files)) {
    for (plot_type in c('density', 'ridgeline', 'histogram')) {
      if (any(grepl(paste0('_', plot_type, '.tsv'), names(pggf$files)))) {
        cli::cli_abort('{.fn pggf_{plot_type}} can only be used as the last element of a {.fn pggfplot}')
      } 
    }
  }
}


# helper to use quoted or unquoted arguments to select columns for ggplot
.pggf_sym_or_str <- function(x) {
  if (is.symbol(x)) {
    return(x)
  } else if (is.null(x)) {
    return(NULL)
  } else if (rlang::is_quosure(x)) {
    x <- rlang::as_name(x)
    return(rlang::quo(.data[[x]]))
  } else {
    return(rlang::quo(.data[[x]]))
  }
}


# helper to test if argument is null
.pggf_is_null <- function(x) {
  return(! is.symbol(substitute(x)) && is.null(x) || (rlang::is_quosure(x) && rlang::quo_is_null(x)))
}


# stats
pggf_stats <- function(pggf, fns) {
  # make sure `pggf` is a `pggf` class object
  .assert_pggf(pggf)
  
  # test if `pggf_stats` has already been called before
  if (! is.null(pggf$stats)) {
    cli::cli_abort(
      c('{.fn pggf_stats} can only be used once per {.fn pggfplot}'),
      i = 'Multiple stats can be combined in {.arg fns}; e.g. "c(!! pggf_stat_fns$correlation, !! pggf_stat_fns$trendline, \'sum\' = ~ sum(x))"')
  }
  
  # test if `fns` is given
  if (missing(fns) || .pggf_is_null(substitute(fns))) {
    cli::cli_abort('argument {.arg fns} is missing, with no default')
  }
  
  # calculate stats
  fns <- rlang::enquo(fns)
  
  stats <- pggf$data |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
    dplyr::summarise(
      dplyr::across(.cols = 1, .fns = !! fns, .names = '{.fn}'),
      .groups = 'drop'
    )
  
  # save stats
  pggf$stats <- readr::format_tsv(stats)
  
  # add to ggplot
  if (! is.null(pggf$plot)) {
    plot_stats <- stats |>
      tidyr::pivot_longer(
        -tidyselect::any_of('facet_id')
      ) |>
      dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
      dplyr::summarise(
        label = paste(paste0(' ', name), round(value, 2), sep = ' = ', collapse = '\n'),
        .groups = 'drop'
      )
    
    if (! is.null(pggf$facets)) {
      plot_stats <- plot_stats |>
        dplyr::inner_join(pggf$facets, by = 'facet_id')
      
      if ('wrap' %in% names(plot_stats)) {
        plot_stats <- plot_stats |> dplyr::rename(!! pggf$facet_wrap := wrap)
      } else {
        if ('col' %in% names(plot_stats)) {plot_stats <- plot_stats |> dplyr::rename(!! pggf$facet_col := col)}
        if ('row' %in% names(plot_stats)) {plot_stats <- plot_stats |> dplyr::rename(!! pggf$facet_row := row)}
      }
    }
    
    pggf$plot <- append(pggf$plot, ggplot2::geom_text(data = plot_stats, mapping = ggplot2::aes(label = label), x = -Inf, y = Inf, hjust = 0, vjust = 1.05))
  }
  
  # add to code
  if (! is.null(pggf$code)) {
    stat_cols <- stats |>
      dplyr::select(-tidyselect::any_of('facet_id')) |>
      names() |>
      paste0(collapse = ',')
    
    pggf$code <- append(pggf$code, paste0('\t\\pggf_stats(stats = {', stat_cols, '})\n'), after = length(pggf$code) - 1)
  }
  
  # call the new `pggf` object
  pggf
}

# pre-defined stat functions
pggf_stat_fns <- list(
  'n' = rlang::quo(c('n' = ~ dplyr::n())),
  'correlation' = rlang::quo(
    c('n' = ~ dplyr::n(), 'pearson' = ~ cor(x, y), 'rsquare' = ~ cor(x, y)^2, 'spearman' = ~ cor(x, y, method = 'spearman'))
  ),
  'trendline' = rlang::quo(
    c('slope' = ~ coef(lm(y ~ x))[2], 'intercept' = ~ coef(lm(y ~ x))[1], 'goodness of fit' = ~ summary(lm(y ~ x))$r.squared)
  )
)


# scatter plot
pggf_scatter <- function(pggf, color = NULL, shuffle = TRUE, split = NULL, extra_cols = NULL, y_error = NULL, x_error = NULL) {
  # make sure `pggf` is a `pggf` class object
  .assert_pggf(pggf, dims = 2)
  
  # prepare `data`
  data <- pggf$data |>
    dplyr::ungroup() |>
    dplyr::select(tidyselect::any_of('facet_id'), x, y, {{color}}, {{split}}, {{extra_cols}}, {{y_error}}, {{x_error}})
  
  # save error data for axis limit calculations
  if (! missing(y_error) && ! .pggf_is_null(substitute(y_error))) {
    y_error_data <- data |>
      dplyr::select({{y_error}})
    
    if (ncol(y_error_data) %in% c(1, 2)) {
      pggf$error[c(3, 4)] <- y_error_data
    } else {
      cli::cli_abort('{.arg y_error} must be of length 1 or 2')
    }
  }
  if (! missing(x_error) && ! .pggf_is_null(substitute(x_error))) {
    x_error_data <- data |>
      dplyr::select({{x_error}})
    
    if (ncol(x_error_data) %in% c(1, 2)) {
      pggf$error[c(1, 2)] <- x_error_data
    } else {
      cli::cli_abort('{.arg x_error} must be of length 1 or 2')
    }
  }
  
  # shuffle data
  if (shuffle) {
    data <- data |>
      dplyr::slice_sample(prop = 1)
  }
  
  # split (optional) and export data
  files <- data |>
    dplyr::arrange(dplyr::pick(tidyselect::any_of('facet_id'))) |>
    dplyr::mutate('.file_basename' = pggf$filename, '.file_suffix' = 'scatter.tsv') |>
    tidyr::unite('file', .file_basename, {{split}}, .file_suffix, remove = FALSE) |>
    dplyr::select(-c(.file_basename, .file_suffix)) |>
    dplyr::nest_by(file) |>
    dplyr::summarise(
      content = readr::format_tsv(data),
      .groups = 'drop'
    ) |>
    tibble::deframe()
  
  # add file(s) to pggf object
  pggf$files <- c(pggf$files, files)
  
  # add to ggplot
  if (! is.null(pggf$plot)) {
    if (! (missing(color) || .pggf_is_null(substitute(color))) && ! (missing(split) || .pggf_is_null(substitute(split)))) {
      plot_color <- .pggf_sym_or_str(substitute(color))
      plot_split <- .pggf_sym_or_str(substitute(split))
      pggf$plot <- append(pggf$plot, ggplot2::geom_point(mapping = ggplot2::aes(color = {{plot_color}}, shape = {{plot_split}})))
    } else if (! missing(color) && ! .pggf_is_null(substitute(color))) {
      plot_color <- .pggf_sym_or_str(substitute(color))
      pggf$plot <- append(pggf$plot, ggplot2::geom_point(mapping = ggplot2::aes(color = {{plot_color}})))
    } else if (! missing(split) && ! .pggf_is_null(substitute(split))) {
      plot_split <- .pggf_sym_or_str(substitute(split))
      pggf$plot <- append(pggf$plot, ggplot2::geom_point(mapping = ggplot2::aes(shape = {{plot_split}})))
    } else {
      pggf$plot <- append(pggf$plot, ggplot2::geom_point())
    }
    if (! missing(y_error) && ! .pggf_is_null(substitute(y_error))) {
      plot_y_err <- names(y_error_data)
      if (length(plot_y_err) == 1) {
        plot_y_err[2] <- plot_y_err
      }
      pggf$plot <- append(pggf$plot, ggplot2::geom_errorbar(ggplot2::aes(ymin = y - .data[[plot_y_err[1]]], ymax = y + .data[[plot_y_err[2]]]), width = 0.1))
    }
    if (! missing(x_error) && ! .pggf_is_null(substitute(x_error))) {
      plot_x_err <- names(x_error_data)
      if (length(plot_x_err) == 1) {
        plot_x_err[2] <- plot_x_err
      }
      pggf$plot <- append(pggf$plot, ggplot2::geom_errorbar(ggplot2::aes(xmin = x - .data[[plot_x_err[1]]], xmax = x + .data[[plot_x_err[2]]]), width = 0.1))
    }
  }
  
  # add to code
  if (! is.null(pggf$code)) {
    if (! missing(y_error) && ! .pggf_is_null(substitute(y_error))) {
      if (length(names(y_error_data)) == 1) {
        code_y_err <- paste0('y error = {', names(y_error_data), '}')
      } else {
        code_y_err <- paste0('y error* = {', names(y_error_data)[1], '}{', names(y_error_data)[2], '}')
      }
    } else {
      code_y_err <- NULL
    }
    if (! missing(x_error) && ! .pggf_is_null(substitute(x_error))) {
      if (length(names(x_error_data)) == 1) {
        code_x_err <- paste0('x error = {', names(x_error_data), '}')
      } else {
        code_x_err <- paste0('x error* = {', names(x_error_data)[1], '}{', names(x_error_data)[2], '}')
      }
    } else {
      code_x_err <- NULL
    }
    
    if (! missing(color) && ! .pggf_is_null(substitute(color))) {
      colors <- data |> dplyr::distinct(dplyr::pick({{color}})) |> dplyr::pull() |> sort() |> paste0(collapse = ',')
      color_str <- if (is.symbol(substitute(color))) {deparse(substitute(color))} else {as.character(color)}
      code_color <- paste0('style from column = ', color_str)
      pggf$code <- append(
        pggf$code,
        paste0('\tscatter styles = {', colors, '},'),
        after = 1
      )
    } else {
      code_color <- NULL
    }
    
    code_extras <- paste0(c(code_x_err, code_y_err, code_color), collapse = ', ')
    
    if (! missing(split) && ! .pggf_is_null(substitute(split))) {
      if (code_extras != '') {
        code_extras <- paste0(', ', code_extras)
      }
      code_split <- paste0('data = ', gsub('_scatter.tsv', '', names(files)), code_extras)
      pggf$code <- append(
        pggf$code,
        c(paste0('\t\\pggf_scatter(', code_split, ')'), '\t'),
        after = length(pggf$code) - 1
      )
    } else {
      pggf$code <- append(
        pggf$code,
        paste0('\t\\pggf_scatter(', code_extras, ')\n'),
        after = length(pggf$code) - 1
      )
    }
    
    if (! is.null(code_color)) {
      pggf$code <- append(pggf$code, paste0('\t\\pggf_annotate(scatter legend)\n'), after = length(pggf$code) - 1)
    }
  }
  
  # call the new `pggf` object
  pggf
}


# box plots
pggf_boxplot <- function(pggf, group = NULL, width = 'fixed', orientation = 'x',
                         extra_cols = NULL, points_below = pggf_defaults$points_below,
                         signif = NULL, type = 'boxplot', trim = pggf_defaults$trim,
                         scale = 'area', samples = pggf_defaults$samples) {
  # make sure `pggf` is a `pggf` class object
  .assert_pggf(pggf, dims = 2)
  
  # check type
  if (! type %in% c('boxplot', 'violin', 'halfviolin')) {
    cli::cli_abort('{.arg type} must be one of "boxplot", "violin", or "halfviolin"')
  }
  
  # get orientation
  if (orientation %in% c('x', 'v', 'vertical')) {
    orientation <- 'x'
    if (! pggf$numeric_y) {
      cli::cli_abort('{.arg y} must be numeric for {.arg orientation = "x"}')
    }
    sample_cols <- c('x', 'xticklabels')
    value_col <- 'y'
  } else if (orientation %in% c('y', 'h', 'horizontal')) {
    orientation <- 'y'
    if (! pggf$numeric_x) {
      cli::cli_abort('{.arg x} must be numeric for {.arg orientation = "y"}')
    }
    sample_cols <- c('y', 'yticklabels')
    value_col <- 'x'
  } else {
    cli::cli_abort(c('{.arg orientation} must be "x" or "y"', i = 'You can also use "v" or "vertical" for "x" and "h" or "horizontal" for "y"'))
  }
  
  # prepare `data`
  box_data <- pggf$data |>
    dplyr::ungroup() |>
    dplyr::mutate(
      dplyr::across({{group}}, as.ordered)
    )
  
  # check if individual points should be drawn
  show_points <- box_data |>
    dplyr::count(dplyr::pick(tidyselect::any_of(c('facet_id', sample_cols)), {{group}})) |>
    dplyr::summarise(any(n < points_below), .groups = 'drop') |>
    dplyr::pull()
  
  # scale box width relative to the group size
  if (! missing(group) && ! .pggf_is_null(substitute(group))) {
    if (width == 'fixed') {
      box_data <- box_data |>
        dplyr::mutate(
          group_n = dplyr::n_distinct(pick({{group}}))
        )
    } else if (width == 'free') {
      box_data <- box_data |>
        dplyr::group_by(dplyr::pick(tidyselect::any_of(c('facet_id', sample_cols)))) |>
        dplyr::mutate(
          group_n = dplyr::n_distinct(pick({{group}}))
        )
    } else {
      cli::cli_abort('{.arg width} must be "fixed" or "free"')
    }
    
    box_data <- box_data |>
      dplyr::mutate(
        box_sep = dplyr::if_else(group_n == 1, 0, 0.1 * sqrt(group_n - 1) / (group_n - 1)),
        box_width = (1 - (group_n - 1) * box_sep) / group_n,
        across({{group}}, ~ (as.numeric(droplevels(.x)) - 0.5 * (group_n + 1)) * (box_width + box_sep), .names = 'box_shift')
      ) |>
      dplyr::ungroup()
  } else {
    box_data <- box_data |>
      dplyr::mutate(
        box_width = 1,
        box_shift = 0
      )
  }
  
  # calculate box plot metrics
  box_metrics <- box_data |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of(c('facet_id', sample_cols)), {{group}}, box_width, box_shift)) |>
    dplyr::summarise(
      n = dplyr::n(),
      med = median(.data[[value_col]], na.rm = TRUE),
      av = mean(.data[[value_col]], na.rm = TRUE),
      lq = quantile(.data[[value_col]], 0.25, na.rm = TRUE),
      uq = quantile(.data[[value_col]], 0.75, na.rm = TRUE),
      iqr = 1.5 * IQR(.data[[value_col]], na.rm = TRUE),
      lw = max(lq - iqr, min(.data[[value_col]][.pggf_outliers(.data[[value_col]]) == 0])),
      uw = min(uq + iqr, max(.data[[value_col]][.pggf_outliers(.data[[value_col]]) == 0])),
      min = min(.data[[value_col]]),
      max = max(.data[[value_col]]),
      outliers = sum(.pggf_outliers(.data[[value_col]]) != 0),
      dplyr::across(
        {{extra_cols}},
        ~ if (length(unique(.x)) == 1) {
          unique(.x)
        } else {
          cli::cli_warn('Non-unique values in {.arg extra_cols}')
          paste0(.x, collapse = ',')
        }
      ),
      .groups = 'drop'
    ) |>
    dplyr::select(-c(iqr)) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::pick(tidyselect::any_of(c('facet_id', orientation)), {{group}})) |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
    dplyr::mutate(
      dplyr::across(tidyselect::any_of(orientation), dplyr::n_distinct, .names = 'n_samples'),
      sample_id = seq_len(dplyr::n())
    ) |>
    dplyr::ungroup()
  
  if (show_points) {
    box_metrics$outliers <- -1
  }
  
  # add significance groups
  if (! is.null(signif)) {
    signif <- as.list(signif)
    
    if (type == 'halfviolin' & is.null(signif$comparisons)) {
      signif$comparisons <- 'sample groups'
    }
    
    if (is.null(pggf$facets)) {
      signif$data <- box_data |>
        dplyr::select(tidyselect::any_of(c(value_col, sample_cols)), {{group}})
    } else {
      signif$data <- box_data |>
        dplyr::select(tidyselect::any_of(c('facet_id', value_col, sample_cols)), {{group}}) |>
        dplyr::inner_join(pggf$facets, by = 'facet_id')
    }
    
    if (pggf[[paste0('numeric_', orientation)]]) {
      signif$data <- signif$data |>
        dplyr::mutate(
          !! sample_cols[2] := .data[[sample_cols[1]]]
        ) |>
        dplyr::select(tidyselect::any_of(c('facet_id', value_col, sample_cols)), {{group}})
    }
    
    signif$data <- signif$data |>
      dplyr::mutate(
        dplyr::across(
          c(tidyselect::any_of(sample_cols[1]), {{group}}),
          as.ordered
        )
      ) |>
      dplyr::arrange(dplyr::pick(tidyselect::any_of(sample_cols), {{group}}))
    
    if (is.null(signif$comparisons) || (is.character(signif$comparisons) && signif$comparisons == 'all')) {# compare all samples with each other; yields compact letter display of significance groups
      signif_type <- 'cld'
      # compute significance groups
      cld <- do.call(
        .pggf_cld,
        signif
      )
      
      # convert numeric samples back
      if (pggf[[paste0('numeric_', orientation)]]) {
        cld <- cld |>
          dplyr::mutate(
            !! sample_cols[1] := as.numeric(.data[[sample_cols[2]]])
          ) |>
          dplyr::select(-tidyselect::all_of(sample_cols[2]))
      }
      
      box_metrics <- box_metrics |>
        dplyr::full_join(
          cld,
          by = names(cld)[names(cld) != 'cld']
        ) 
    } else if (is.character(signif$comparisons) && signif$comparisons == 'one-sample tests') {# one-sample tests; yield signifcance labels
      signif_type <- 'label'
      
      signif_labels <- do.call(
        .pggf_onesample_tests,
        signif
      )
      
      # convert numeric samples back
      if (pggf[[paste0('numeric_', orientation)]]) {
        signif_labels <- signif_labels |>
          dplyr::mutate(
            !! sample_cols[1] := as.numeric(.data[[sample_cols[2]]])
          ) |>
          dplyr::select(-tidyselect::all_of(sample_cols[2]))
      }
      
      box_metrics <- box_metrics |>
        dplyr::full_join(
          signif_labels,
          by = names(signif_labels)[names(signif_labels) != 'p_value']
        ) 
    } else {# only compare specified samples; yield significance brackets
      signif_type <- 'brackets'
      # compute significance brackets
      signif_brackets <- do.call(
        .pggf_signif,
        signif
      )
      
      # unify max and range across facets
      if (! is.null(pggf$facets) && pggf$scales %in% c('free', paste0('free_', value_col))) {
        if (orientation == 'x') {
          signif_brackets <- signif_brackets |>
            dplyr::group_by(dplyr::pick(tidyselect::any_of('row'))) |>
            dplyr::mutate(
              min_y = min(min_y),
              max_y = max(max_y)
            )
        } else {
          signif_brackets <- signif_brackets |>
            dplyr::group_by(dplyr::pick(tidyselect::any_of('col'))) |>
            dplyr::mutate(
              min_x = min(min_x),
              max_x = max(max_x)
            )
        }
        
        signif_brackets <- signif_brackets |>
          dplyr::ungroup() |>
          dplyr::select(-tidyselect::any_of(c('col', 'row', 'wrap')))
      } else {
        signif_brackets <- signif_brackets |>
          dplyr::mutate(
            !! paste0('min_', value_col) := min(.data[[paste0('min_', value_col)]]),
            !! paste0('max_', value_col) := max(.data[[paste0('max_', value_col)]])
          )
      }
      
      signif_brackets <- signif_brackets |>
        dplyr::mutate(
          range = .data[[paste0('max_', value_col)]] - .data[[paste0('min_', value_col)]]
        ) |>
        dplyr::select(-tidyselect::all_of(paste0('min_', value_col))) |>
        dplyr::relocate(range, .after = tidyselect::all_of(paste0('max_', value_col)))
      
      # export data
      signif_file <- signif_brackets |>
        dplyr::mutate(
          file = paste0(pggf$filename, '_signif.tsv')
        ) |>
        dplyr::nest_by(file) |>
        dplyr::summarise(
          content = readr::format_tsv(data),
          .groups = 'drop'
        ) |>
        tibble::deframe()
      
      # add file to pggf object
      pggf$files <- c(pggf$files, signif_file)
    }
  }
  
  # export data
  box_file <- box_metrics |>
    dplyr::mutate(
      file = paste0(pggf$filename, '_boxplot.tsv')
    ) |>
    dplyr::nest_by(file) |>
    dplyr::summarise(
      content = readr::format_tsv(data),
      .groups = 'drop'
    ) |>
    tibble::deframe()
  
  # add file to pggf object
  pggf$files <- c(pggf$files, box_file)
  
  # create files for violin or for points/outliers
  if (type %in% c('violin', 'halfviolin')) {
    # prepare `data` for violin plots
    sample_ids <- box_metrics |>
      dplyr::select(tidyselect::any_of(c('sample_id', 'facet_id', orientation)), {{group}})
    
    violin_data <- pggf$data |>
      dplyr::mutate(
        dplyr::across({{group}}, as.ordered)
      ) |>
      dplyr::inner_join(sample_ids, by = names(sample_ids)[-1])
    
    # calculate density
    violin_data <- violin_data |>
      dplyr::ungroup() |>
      dplyr::nest_by(dplyr::pick(tidyselect::any_of(c('facet_id', 'sample_id')), {{group}})) |>
      dplyr::reframe(
        density = dplyr::if_else(
          trim,
          list(density(data[[value_col]], from = min(data[[value_col]]), to = max(data[[value_col]]), n = {{samples}})),
          list(density(data[[value_col]], n = {{samples}}))
        ),
        x = density[[value_col]],
        y = density[[orientation]],
        n = length(data[[value_col]])
      ) |>
      dplyr::ungroup() |>
      dplyr::select(-density)
    
    # scale width
    if (scale == 'area') {
      violin_data <- violin_data |>
        dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
        dplyr::mutate(
          {{orientation}} := 0.5 * .data[[orientation]] / max(.data[[orientation]])
        ) |>
        dplyr::ungroup()
    } else if (scale == 'width') {
      violin_data <- violin_data |>
        dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id')), sample_id, {{group}}) |>
        dplyr::mutate(
          {{orientation}} := 0.5 * .data[[orientation]] / max(.data[[orientation]])
        ) |>
        dplyr::ungroup()
    } else if (scale == 'count') {
      violin_data <- violin_data |>
        dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
        dplyr::mutate(
          {{orientation}} := 0.5 * (.data[[orientation]] * n) / (max(.data[[orientation]] * n))
        ) |>
        dplyr::ungroup()
    } else {
      cli::cli_abort('{.arg scale} must be one of "area", "width", or "count"')
    }
    
    violin_data <- violin_data  |>
      dplyr::select(-n) |>
      dplyr::arrange(dplyr::pick(tidyselect::any_of('facet_id')), sample_id, {{group}})
    
    # export data
    violin_file <- violin_data |>
      dplyr::mutate(
        file = paste0(pggf$filename, '_violin.tsv')
      ) |>
      dplyr::nest_by(file) |>
      dplyr::summarise(
        content = readr::format_tsv(data),
        .groups = 'drop'
      ) |>
      tibble::deframe()
    
    # add file to pggf object
    pggf$files <- c(pggf$files, violin_file)
  } else if (show_points || any(box_metrics$outliers > 0)) {
    # create file for points or outliers
    outliers <- box_data |>
      dplyr::ungroup() |>
      dplyr::select(tidyselect::any_of(c('facet_id', sample_cols, value_col)), box_shift, {{group}}, {{extra_cols}})
    
    if (! show_points) {
      outliers <- outliers |>
        dplyr::group_by(dplyr::pick(tidyselect::any_of(c('facet_id', sample_cols)), {{group}})) |>
        dplyr::filter(.pggf_outliers(.data[[value_col]]) != 0) |>
        dplyr::ungroup()
    }
    
    if (missing(group) || .pggf_is_null(substitute(group))) {
      if ('facet_id' %in% names(box_metrics)) {
        outliers <- outliers |>
          dplyr::inner_join(
            box_metrics |> dplyr::select(tidyselect::any_of(c('facet_id', orientation)), sample_id),
            dplyr::join_by('facet_id', {{orientation}})
          )
      } else {
        outliers <- outliers |>
          dplyr::inner_join(
            box_metrics |> dplyr::select(tidyselect::any_of(c('facet_id', orientation)), sample_id),
            dplyr::join_by({{orientation}})
          ) 
      }
    } else {
      if ('facet_id' %in% names(box_metrics)) {
        outliers <- outliers |>
          dplyr::inner_join(
            box_metrics |> dplyr::select(tidyselect::any_of(c('facet_id', orientation)), {{group}}, sample_id),
            dplyr::join_by('facet_id', {{orientation}}, {{group}})
          )
      } else {
        outliers <- outliers |>
          dplyr::inner_join(
            box_metrics |> dplyr::select(tidyselect::any_of(c('facet_id', orientation)), {{group}}, sample_id),
            dplyr::join_by({{orientation}}, {{group}})
          ) 
      }
    }
    
    outliers <- outliers |>
      dplyr::mutate(
        dplyr::across(tidyselect::ends_with('ticklabels'), ~ dplyr::if_else(grepl(' ', .x, fixed = TRUE), paste0('{', .x, '}'), .x)),
        file = paste0(pggf$filename, '_outlier.tsv')
      ) |>
      dplyr::nest_by(file) |>
      dplyr::summarise(
        content = readr::format_tsv(data),
        .groups = 'drop'
      ) |>
      tibble::deframe()
    
    pggf$files <- c(pggf$files, outliers)
  }
  
  # add to ggplot
  if (! is.null(pggf$plot)) {
    if (type %in% c('violin', 'halfviolin')) {
      plot_width <- 0.1
    } else {
      plot_width <- 0.75
    }
    
    if (missing(group) || .pggf_is_null(substitute(group))) {
      if (type %in% c('violin', 'halfviolin')) {
        pggf$plot <- append(
          pggf$plot,
          ggplot2::geom_violin(
            mapping = ggplot2::aes(group = .data[[orientation]]),
            orientation = orientation,
            scale = scale,
            trim = trim
          )
        )
      }
      pggf$plot <- append(
        pggf$plot,
        c(
          ggplot2::geom_boxplot(
            mapping = ggplot2::aes(group = .data[[orientation]]),
            orientation = orientation,
            outlier.shape = ifelse(show_points || type %in% c('violin', 'halfviolin'), NA, 19),
            width = plot_width
          ),
          ggplot2::stat_summary(
            fun.data = .pggf_sample_size,
            geom = 'text',
            orientation = orientation,
            hjust = dplyr::if_else(orientation == 'x', 0.5, -0.5),
            vjust = dplyr::if_else(orientation == 'x', -1, 0.5)
          )
        )
      )
    } else {
      plot_group <- .pggf_sym_or_str(substitute(group))
      if (type %in% c('violin', 'halfviolin')) {
        pggf$plot <- append(
          pggf$plot,
          ggplot2::geom_violin(
            mapping = ggplot2::aes(group = interaction(.data[[orientation]], {{plot_group}}), fill = as.ordered({{plot_group}})),
            orientation = orientation,
            scale = scale,
            trim = trim
          )
        )
      }
      pggf$plot <- append(
        pggf$plot,
        c(
          ggplot2::geom_boxplot(
            mapping = ggplot2::aes(group = interaction(.data[[orientation]], {{plot_group}}), fill = as.ordered({{plot_group}})),
            orientation = orientation,
            outlier.shape = ifelse(show_points || type %in% c('violin', 'halfviolin'), NA, 19),
            width = plot_width,
            position = ggplot2::position_dodge(0.9)
          ),
          ggplot2::stat_summary(
            fun.data = .pggf_sample_size,
            ggplot2::aes(color = as.ordered({{plot_group}})),
            geom = 'text',
            orientation = orientation,
            hjust = dplyr::if_else(orientation == 'x', 0.5, -0.5),
            vjust = dplyr::if_else(orientation == 'x', -1, 0.5),
            position = ggplot2::position_dodge(0.9)
          )
        )
      )
    }
    
    if (show_points) {
      if (missing(group) || .pggf_is_null(substitute(group))) {
        pggf$plot <- append(pggf$plot, ggplot2::geom_jitter(width = 0.3))
      } else {
        pggf$plot <- append(
          pggf$plot,
          ggplot2::geom_point(
            mapping = ggplot2::aes(fill = as.ordered({{plot_group}})),
            position = ggplot2::position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9)
          )
        )
      }
    }
    
    if (! is.null(signif)) {
      if (signif_type == 'cld') {
        plot_cld <- pggf$data |>
          dplyr::group_by(
            dplyr::pick(
              tidyselect::any_of(c('facet_id', sample_cols)),
              {{group}},
              !! pggf$facet_col,
              !! pggf$facet_row,
              !! pggf$facet_wrap
            )
          ) |>
          dplyr::summarise(
            dplyr::across(tidyselect::all_of(value_col), max),
            .groups = 'drop'
          ) |>
          dplyr::mutate(
            dplyr::across({{group}}, as.ordered)
          )
        
        if (missing(group) || .pggf_is_null(substitute(group))) {
          plot_cld <- plot_cld |>
            dplyr::inner_join(cld, by = names(cld)[names(cld) != 'cld'])
          
          pggf$plot <- append(
            pggf$plot,
            ggplot2::geom_text(
              data = plot_cld,
              mapping = ggplot2::aes(label = cld),
              hjust = dplyr::if_else(orientation == 'x', 0.5, -1),
              vjust = dplyr::if_else(orientation == 'x', -1, 0.5)
            )
          )
        } else {
          plot_cld <- cld |>
            dplyr::mutate(
              {{group}} := ordered({{group}}, levels = plot_cld |> dplyr::pull({{group}}) |> levels())
            ) |>
            dplyr::inner_join(plot_cld, by = names(cld)[names(cld) != 'cld'])
          
          pggf$plot <- append(
            pggf$plot,
            ggplot2::geom_text(
              data = plot_cld,
              mapping = ggplot2::aes(label = cld, group = interaction(.data[[orientation]], {{plot_group}})),
              position = ggplot2::position_dodge(0.9),
              hjust = dplyr::if_else(orientation == 'x', 0.5, -1),
              vjust = dplyr::if_else(orientation == 'x', -1, 0.5)
            )
          )
        }
      } else if (signif_type == 'label') {
        plot_signif <- pggf$data |>
          dplyr::group_by(
            dplyr::pick(
              tidyselect::any_of(c('facet_id', sample_cols)),
              {{group}},
              !! pggf$facet_col,
              !! pggf$facet_row,
              !! pggf$facet_wrap
            )
          ) |>
          dplyr::summarise(
            dplyr::across(tidyselect::all_of(value_col), max),
            .groups = 'drop'
          ) |>
          dplyr::mutate(
            dplyr::across({{group}}, as.ordered)
          ) |>
          dplyr::inner_join(signif_labels, by = names(signif_labels)[names(signif_labels) != 'p_value'])
        
        if (missing(group) || .pggf_is_null(substitute(group))) {
          pggf$plot <- append(
            pggf$plot,
            ggplot2::geom_text(
              data = plot_signif,
              mapping = ggplot2::aes(label = p_value),
              hjust = dplyr::if_else(orientation == 'x', 0.5, -1),
              vjust = dplyr::if_else(orientation == 'x', -1, 0.5)
            )
          )
        } else {
          pggf$plot <- append(
            pggf$plot,
            ggplot2::geom_text(
              data = plot_signif,
              mapping = ggplot2::aes(label = p_value, group = interaction(.data[[orientation]], {{plot_group}})),
              position = ggplot2::position_dodge(0.9),
              hjust = dplyr::if_else(orientation == 'x', 0.5, -1),
              vjust = dplyr::if_else(orientation == 'x', -1, 0.5)
            )
          )
        }
      }
    }
  }
  
  # add to code
  if (! is.null(pggf$code)) {
    if (! is.null(signif)) {
      if (signif_type == 'cld') {
        code_enlarge <- paste0('\tenlarge ', ifelse(orientation == 'x', 'y', 'x'), ' limits = {lower = 0.1, upper = 0.1}, % increase axis limits to fit sample sizes and cld')
        code_cld <- 'cld = {anchor = south, text depth = 0pt}{max}'
        code_signif <- NULL
      } else if (signif_type == 'label') {
        code_enlarge <- paste0('\tenlarge ', ifelse(orientation == 'x', 'y', 'x'), ' limits = {lower = 0.1, upper = 0.1}, % increase axis limits to fit sample sizes and sigificance label')
        code_cld <- 'p value = {anchor = south, text depth = 0pt}{max}'
        code_signif <- NULL
      } else {
        code_enlarge <- paste0('\tenlarge ', ifelse(orientation == 'x', 'y', 'x'), ' limits = {lower = 0.1, upper = 0.2}, % increase axis limits to fit sample sizes and significance brackets')
        code_cld <- NULL
        code_signif <- '\t\\pggf_signif()\n'
      }
    } else {
      code_enlarge <- paste0('\tenlarge ', ifelse(orientation == 'x', 'y', 'x'), ' limits = {lower = 0.1}, % increase axis limits to fit sample sizes')
      code_cld <- NULL
      code_signif <- NULL
    }
    
    if (type == 'halfviolin') {
      code_extra <- paste0(c('half violin', code_cld), collapse = ', ')
    } else {
      code_extra <- code_cld
    }
    
    pggf$code <- append(
      pggf$code,
      code_enlarge,
      after = 1
    )
    pggf$code <- append(pggf$code, paste0('\t\\pggf_', ifelse(type == 'halfviolin', 'violin', type), '(', code_extra, ')\n'), after = length(pggf$code) - 1)
    
    if (! is.null(code_signif)) {
      pggf$code <- append(pggf$code, code_signif, after = length(pggf$code) - 1)
    }
  }
  
  # call the new `pggf` object
  pggf
}


# find outliers
.pggf_outliers <- function(data) {
  dplyr::case_when(
    data > quantile(data, 0.75, na.rm = TRUE) + 1.5 * IQR(data, na.rm = TRUE) ~ 1,
    data < quantile(data, 0.25, na.rm = TRUE) - 1.5 * IQR(data, na.rm = TRUE) ~ -1,
    TRUE ~ 0
  )
}

# display sample size on ggplot
.pggf_sample_size <- function(x, y = -Inf, size = 8/ggplot2::.pt){
  return(c(y = y, size = size, label = length(x)))
}


# violin plots
pggf_violin <- function(pggf, half = NULL, group = NULL, orientation = 'x', extra_cols = NULL,
                        width = 'fixed', trim = pggf_defaults$trim, scale = 'area',
                        samples = pggf_defaults$samples, signif = NULL) {
  # check arguments
  if (! is.logical(trim) || length(trim) != 1) {
    cli::cli_abort('{.arg trim} must be a single logical value.')
  }
  if (! is.numeric(samples) || length(samples) != 1) {
    cli::cli_abort('{.arg samples} must be a single number.')
  }
  if (! missing(half) && ! .pggf_is_null(substitute(half))) {
    if (! missing(group) && ! .pggf_is_null(substitute(group))) {
      cli::cli_abort('{.arg half} and {.arg group} are mutually exclusive. Only one can be used at a time.')
    }
    
    half_levels <- pggf$data |>
      dplyr::distinct({{half}}) |>
      nrow()
    
    if (half_levels == 2) {
      plottype = 'halfviolin'
      group = rlang::enquo(half)
    } else {
      cli::cli_abort('{.arg half} must specify a column of the data with exactly two distinct values.')
    }
  } else {
    plottype = 'violin'
  }
  
  # call pggf_boxplot
  pggf <- rlang::inject(
    pggf_boxplot(
      pggf,
      group = {{group}},
      orientation = orientation,
      extra_cols = {{extra_cols}},
      width = width,
      type = plottype,
      points_below = 0,
      trim = trim,
      scale = scale,
      samples = samples,
      signif = signif
    )
  )
  
  # call the new `pggf` object
  pggf
}


# significance groups in compact letter display
.pggf_cld <- function(data, test, threshold = 0.05, p_adjust = NULL, order_by = 'sample', 
                     reverse_order = FALSE, save_test_results = NULL, ...) {
  # check arguments
  if (! is.null(p_adjust)) {
    if (! p_adjust %in% p.adjust.methods) {
      cli::cli_abort('{.arg p_adjust} must be one of: {p.adjust.methods}')
    }
    if (test == 'Tukey') {
      cli::cli_warn('Selecting a p value adjustment method is not possible for the Tukey test. Ignoring {.arg p_adjust}.')
    }
  }
  
  # get column names and prepare data
  facet_cols <- names(data)[names(data) %in% c('facet_id', 'col', 'row', 'wrap')]
  col_names <- names(data)[! names(data) %in% c('facet_id', 'col', 'row', 'wrap')]
  value <- col_names[1]
  sample_order <- col_names[2]
  group <- col_names[-c(1, 2)]
  
  test_data <- data |>
    dplyr::arrange(dplyr::pick(tidyselect::all_of(c(sample_order, group)))) |>
    tidyr::unite(
      'groups',
      tidyselect::all_of(group),
      sep = ':'
    )
  
  group_names <- unique(test_data$groups)
  group_id <- setNames(seq_along(group_names), group_names)
  
  test_data <- test_data |>
    dplyr::mutate(
      groups = as.ordered(group_id[groups])
    ) |>
    dplyr::group_by(dplyr::pick(tidyselect::all_of(facet_cols)))
  
  # run test
  if (test %in% c('t-test', 'Wilcox')) {
    if (test == 't-test') {
      test_results <- test_data |>
        dplyr::summarise(
          test_result = list(pairwise.t.test(.data[[value]], groups, p.adjust.method = p_adjust, ...)),
          .groups = 'keep'
        )
    } else if (test == 'Wilcox') {
      test_results <- test_data |>
        dplyr::summarise(
          test_result = list(pairwise.wilcox.test(.data[[value]], groups, p.adjust.method = p_adjust, ...)),
          .groups = 'keep'
        )
    }
    
    pvals <- test_results |>
      dplyr::mutate(
        p_value = list(
          unlist(test_result, recursive = FALSE)$p.value |>
            tibble::as_tibble(rownames = 'V1') |>
            tidyr::pivot_longer(
              -V1,
              names_to = 'V2',
              values_to = 'p_value',
              values_drop_na = TRUE
            ) |>
            dplyr::arrange(V1, V2) |>
            tidyr::unite(
              'name',
              V1,
              V2,
              sep = '-'
            ) |>
            tibble::deframe()
        )
      )
  } else if (test == 'Tukey') {
    test_results <- test_data |>
      dplyr::summarise(
        anova = list(aov(as.formula(paste0(value, ' ~ groups')))),
        test_result = list(
          anova[[1]] |>
            TukeyHSD(...)
        ),
        anova = summary(anova[[1]])[[1]][1, 'Pr(>F)'],
        .groups = 'keep'
      )
    
    pvals <- test_results |>
      dplyr::summarise(
        p_value = list(
          unlist(test_result, recursive = FALSE)[['groups']][, 'p adj']
        ),
        .groups = 'keep'
      )
  } else {
    cli::cli_abort('Unknown {.arg test} method {.str {test}}. Should be one of "t-test", "Wilcox", or "Tukey"')
  }
  
  # order p-values
  if (order_by == 'sample') {
    group_order <- data |>
      dplyr::ungroup() |>
      dplyr::distinct(dplyr::pick(tidyselect::any_of(c('facet_id', sample_order, group)))) |>
      dplyr::arrange(dplyr::pick(tidyselect::any_of(c('facet_id', sample_order, group[-1])))) |>
      tidyr::unite(
        'groups',
        tidyselect::all_of(group),
        sep = ':'
      ) |>
      dplyr::mutate(
        groups = as.ordered(group_id[groups])
      ) |>
      dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
      dplyr::summarise(
        group_order = list(groups),
        .groups = 'drop'
      )
  } else if (order_by %in% c('median', 'mean')) {
    group_order <- test_data |>
      dplyr::group_by(groups, .add = TRUE) |>
      dplyr::summarise(
        dplyr::across(tidyselect::all_of(value), ~ eval(parse(text = order_by))(.x)),
        .groups = 'drop'
      ) |>
      dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
      dplyr::arrange(dplyr::pick(tidyselect::all_of(value)), .by_group = TRUE) |>
      dplyr::summarise(
        group_order = list(groups),
        .groups = 'drop'
      )
  } else {
    cli::cli_abort('{.arg order_by} must be one of "sample", "median", or "mean"')
  }
  
  if ('facet_id' %in% names(data)) {
    pvals <- pvals |>
      dplyr::inner_join(group_order, by = 'facet_id')
  } else {
    pvals <- pvals |>
      dplyr::mutate(
        group_order = group_order
      )
  }
  
  if (reverse_order) {
    pvals <- pvals |>
      dplyr::mutate(
        p_value = list(
          .pggf_cld_order(unlist(p_value), rev(unlist(group_order)))
        )
      )
  } else {
    pvals <- pvals |>
      dplyr::mutate(
        p_value = list(
          .pggf_cld_order(unlist(p_value), unlist(group_order))
        )
      )
  }
  
  # convert p-values into significance groups
  cld <- pvals |>
    dplyr::reframe(
      p_value |>
        unlist() |>
        multcompView::multcompLetters(threshold = threshold) |>
        _$Letters |>
        tibble::enframe(name = 'groups', value = 'cld')
    ) |>
    dplyr::mutate(
      groups = group_names[as.integer(groups)]
    ) |>
    dplyr::select(tidyselect::any_of('facet_id'), groups, cld) |>
    tidyr::separate_wider_delim(
      groups,
      delim = ':',
      names = group
    )
  
  # save test results
  if (! is.null(save_test_results)) {
    # extract p values and additional information
    if (test %in% c('t-test', 'Wilcox')) {
      test_results <- test_results |>
        dplyr::group_by(dplyr::pick(tidyselect::any_of(c('facet_id', 'col', 'row', 'wrap')))) |>
        dplyr::reframe(
          test = test_result[[1]]$method,
          p_adjust_method = test_result[[1]]$p.adjust.method,
          test_result[[1]]$p.value |>
            tibble::as_tibble(rownames = 'sample2') |>
            tidyr::pivot_longer(
              -sample2,
              names_to = 'sample1',
              values_to = 'p_value',
              values_drop_na = TRUE
            )
        )
    } else if (test == 'Tukey') {
      test_results <- test_results |>
        dplyr::group_by(dplyr::pick(tidyselect::any_of(c('facet_id', 'col', 'row', 'wrap')))) |>
        dplyr::reframe(
          conf.level = test_result[[1]] |> attr('conf.level'),
          test_result[[1]]$groups |>
            tibble::as_tibble(rownames = 'comparison') |>
            dplyr::select(comparison, p_value = `p adj`)
        ) |>
        tidyr::separate_wider_delim(
          comparison,
          delim = '-',
          names = c('sample2', 'sample1')
        ) |>
        dplyr::mutate(
          test = paste0('Tukey multiple comparisons of means (', conf.level * 100,'% family-wise confidence level)'),
        ) |>
        dplyr::bind_rows(
          test_results |>
            dplyr::mutate(
              test = 'ANOVA',
              p_value = anova
            )
        )
    }
    
    # convert samples names back and select columns
    test_results <- test_results |>
      dplyr::arrange(dplyr::pick(tidyselect::any_of(c('facet_id', 'col', 'row', 'wrap')), test, sample1, sample2)) |>
      dplyr::mutate(
        dplyr::across(c(sample1, sample2), ~ group_names[as.integer(.x)])
      ) |>
      tidyr::unite(
        'comparison',
        sample1,
        sample2,
        sep = ' vs. ',
        na.rm = TRUE
      ) |>
      dplyr::select(tidyselect::any_of(c('facet_id', 'col', 'row', 'wrap')), test, comparison, p_value, tidyselect::any_of('p_adjust_method'))
    
    # save to `.pggf_test_results` environment
    assign(save_test_results, test_results |> dplyr::ungroup(), envir = .pggf_test_results)
  }
  
  # return signifcance groups
  return(cld)
}

.pggf_test_results <- new.env()


# function to get saved test results
pggf_get_test_results <- function(names = NULL) {
  if (missing(names) || is.null(names)) {
    names <- ls(envir = .pggf_test_results)
  }
  
  results <- names |>
    sapply(
      function(x) {
        get(x, envir = .pggf_test_results) |>
          dplyr::mutate(
            dplyr::across(tidyselect::where(is.factor), as.character)
          )
      },
      simplify = FALSE,
      USE.NAMES = TRUE
    ) |>
    dplyr::bind_rows(.id = 'name')
  
  return(results)
}


# helper function to order significance groups
.pggf_cld_order <- function(pvals, ord) {
  pvals <- pvals |>
    tibble::enframe(name = 'comparison', value = 'p_value') |>
    tidyr::separate_wider_delim(
      comparison,
      delim = '-',
      names = c('V1', 'V2')
    ) |>
    dplyr::mutate(
      dplyr::across(c(V1, V2), ~ ordered(.x, levels = ord)),
      tmp = dplyr::if_else(V1 > V2, V1, V2),
      V1 = dplyr::if_else(V1 > V2, V2, V1),
      V2 = tmp
    ) |>
    dplyr::select(-tmp) |>
    dplyr::arrange(V1, V2) |>
    tidyr::unite(
      'comparison',
      V1,
      V2,
      sep = '-'
    ) |>
    tibble::deframe()
  
  return(pvals)
}


# significance brackets
.pggf_signif <- function(data, test, comparisons, threshold = 0.05, p_adjust = NULL,
                         save_test_results = NULL, ...) {
  # check arguments
  if (is.character(comparisons) && comparisons == 'sample groups') {
    comp_type <- 'sample groups'
  } else if (is.list(comparisons)) {
    comp_type <- 'list'
    if (is.data.frame(comparisons)) {
      if (! 'sample1' %in% colnames(comparisons) || ! 'sample2' %in% colnames(comparisons)) {
        cli::cli_abort('The {.arg comparisons} dataframe must contain at least the columns "sample1" and "sample2".')
      }
      comp_type <- 'df'
    } else if (any(sapply(comparisons, length) != 2)) {
      cli::cli_abort('Every list element of {.arg comparisons} must be of length 2.')
    }
  }  else {
    cli::cli_abort(
      c(
        '{.arg comparisons} must be a list or a dataframe, not a {.type {comparisons}}',
        'i' = 'You can also use the string "sample groups" to compare all groups of a given sample, the string "all" to compare all samples with each other, or the string "one-sample tests" to compare all samples to a Null distribution.'
        )
      )
  }
  if (is.null(p_adjust)) {
    p_adjust <- p.adjust.methods[1]
  } else if (! p_adjust %in% p.adjust.methods) {
    cli::cli_abort('{.arg p_adjust} must be one of: {p.adjust.methods}')
  }
  
  # set up test function
  if (test == 't-test') {
    test_fn = t.test
  } else if (test == 'Wilcox') {
    test_fn = wilcox.test
  } else {
    cli::cli_abort('Unknown {.arg test} method {.str {test}}. Should be "t-test" or "Wilcox"')
  }
  
  # get column names and prepare data
  facet_cols <- names(data)[names(data) %in% c('facet_id', 'col', 'row', 'wrap')]
  col_names <- names(data)[! names(data) %in% c('facet_id', 'col', 'row', 'wrap')]
  value <- col_names[1]
  sample_order <- col_names[2]
  group <- col_names[-c(1, 2)]
  
  signif_data <- data |>
    dplyr::arrange(dplyr::pick(tidyselect::all_of(sample_order))) |>
    tidyr::unite(
      'groups',
      tidyselect::all_of(group),
      sep = ':'
    )
  
  group_names <- unique(signif_data$groups)
  
  if ('facet_id' %in% colnames(signif_data)) {
    facet_group_names <- signif_data |>
      dplyr::distinct(facet_id, groups) |>
      dplyr::group_by(facet_id) |>
      dplyr::summarise(
        groups = list(groups),
        .groups = 'drop'
      ) |>
      tibble::deframe()
  }
  
  # convert numeric comparisons to sample names
  if (comp_type == 'sample groups') {
    comparisons <- signif_data |>
      dplyr::group_by(dplyr::pick(tidyselect::any_of(c('facet_id', sample_order)))) |>
      dplyr::filter(dplyr::n_distinct(groups) > 1) |>
      dplyr::reframe(
        combn(unique(groups), 2) |>
          t() |>
          tibble::as_tibble(.name_repair = ~ paste0('sample', seq_along(.x)))
      ) |>
      dplyr::select(tidyselect::any_of('facet_id'), sample1, sample2)
    
    if (nrow(comparisons) == 0) {
      cli::cli_abort('{.arg comparisons = "sample groups"} is invalid because no groups are defined.')
    }
  } else if (comp_type == 'list') {
    comparisons <- lapply(comparisons, .pggf_id2sample, group_names) |>
      tibble::tibble() |>
      dplyr::rename('sample' = 1) |>
      tidyr::unnest_wider(sample, names_sep = '')
  } else if (comp_type == 'df') {
    comparisons <- comparisons |>
      dplyr::mutate(
        dplyr::across(c(sample1, sample2), ~ .pggf_id2sample(.x, group_names)),
      )
  }
  
  # add comparisons to data
  if (any(facet_cols %in% colnames(comparisons))) {
    comp_facet_cols <- facet_cols[facet_cols %in% colnames(comparisons)]
    
    signif_data <- signif_data |>
      dplyr::inner_join(comparisons, by = comp_facet_cols, relationship = 'many-to-many')
  } else {
    signif_data <- signif_data |>
      dplyr::cross_join(comparisons)
  }
  
  # check sample names
  bad_samples <- unique(c(signif_data$sample1, signif_data$sample2))
  bad_samples <- bad_samples[! bad_samples %in% unique(signif_data$groups)]
  
  if (length(bad_samples) > 0) {
    cli::cli_abort('Unknown sample names in {.arg comparisons}: {bad_samples}')
  }
  
  if (any(signif_data$sample1 == signif_data$sample2)) {
    cli::cli_abort('{.arg comparisons} contains at least one comparison with identical samples.')
  }
  
  # calculate p values
  signif_data <- signif_data |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
    dplyr::filter(sample1 %in% groups & sample2 %in% groups) |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of(facet_cols)), sample1, sample2) |>
    dplyr::summarise(
      test = list(test_fn(.data[[value]][groups == sample1], .data[[value]][groups == sample2], ...)),
      p_value = test[[1]]$p.value,
      test = test[[1]]$method,
      !! paste0('min_', value) := min(.data[[value]]),
      !! paste0('max_', value) := max(.data[[value]]),
      .groups = 'drop'
    ) |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of(facet_cols))) |>
    dplyr::mutate(
      p_value = p.adjust(p_value, method = p_adjust),
      p_adjust_method = p_adjust,
      dplyr::across(c(sample1, sample2), ~ ordered(.x, levels = group_names)),
      .pggf.tmp = sample2,
      sample2 = dplyr::if_else(as.numeric(sample1) < as.numeric(sample2), sample2, sample1),
      sample1 = dplyr::if_else(as.numeric(sample1) < as.numeric(sample2), sample1, .pggf.tmp),
      .pggf.sample_dist = as.numeric(droplevels(sample2)) - as.numeric(droplevels(sample1))
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::pick(tidyselect::any_of(facet_cols), .pggf.sample_dist, sample1))
  
  # save test results
  if (! is.null(save_test_results)) {
    # order, combine samples and select columns
    test_results <- signif_data |>
      dplyr::arrange(dplyr::pick(tidyselect::any_of(facet_cols), sample1, sample2)) |>
      tidyr::unite(
        'comparison',
        sample1,
        sample2,
        sep = ' vs. ',
        na.rm = TRUE
      ) |>
      dplyr::select(tidyselect::any_of(facet_cols), test, comparison, p_value, tidyselect::any_of('p_adjust_method'))
    
    # save to `.pggf_test_results` environment
    assign(save_test_results, test_results |> dplyr::ungroup(), envir = .pggf_test_results)
  }
  
  # select final columns
  signif_data <- signif_data |>
    dplyr::select(tidyselect::any_of(c(facet_cols, 'sample1', 'sample2', paste0(c('max_', 'min_'), value), 'p_value')))
  
  # convert samples to numeric
  if ('facet_id' %in% colnames(signif_data)) {
    signif_data <- signif_data |>
      dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
      dplyr::mutate(
        dplyr::across(c(sample1, sample2), ~ as.numeric(ordered(.x, levels = facet_group_names[[as.character(unique(facet_id))]])))
      ) |>
      dplyr::ungroup()
  } else {
    signif_data <- signif_data |>
      dplyr::mutate(
        dplyr::across(c(sample1, sample2), as.numeric)
      )
  }
  
  return(signif_data)
}


# helper function to convert sample ids to sample names
.pggf_id2sample <- function(ids, sample_names) {
  if (is.numeric(ids)) {
    samples <- sample_names[ids]
    
    if (any(is.na(samples))) {
      cli::cli_abort('Invalid sample indices in {.arg comparisons}: {ids[is.na(samples)]}')
    }
  } else {
    samples <- ids
  }
  
  return(samples) 
}


# one-sample significance tests
.pggf_onesample_tests <- function(data, test, threshold = 0.05, p_adjust = NULL, save_test_results = NULL, ...) {
  # check arguments
  if (is.null(p_adjust)) {
    p_adjust <- p.adjust.methods[1]
  } else if (! p_adjust %in% p.adjust.methods) {
    cli::cli_abort('{.arg p_adjust} must be one of: {p.adjust.methods}')
  }  
  
  # set up test function
  if (test == 't-test') {
    test_fn = t.test
  } else if (test == 'Wilcox') {
    test_fn = wilcox.test
  } else {
    cli::cli_abort('Unknown {.arg test} method {.str {test}}. Should be "t-test" or "Wilcox"')
  }  
  
  # get column names and prepare data
  facet_cols <- names(data)[names(data) %in% c('facet_id', 'col', 'row', 'wrap')]
  col_names <- names(data)[! names(data) %in% c('facet_id', 'col', 'row', 'wrap')]
  value <- col_names[1]
  sample_order <- col_names[2]
  group <- col_names[-c(1, 2)]
  
  signif_data <- data |>
    tidyr::unite(
      'groups',
      tidyselect::all_of(group),
      sep = ':'
    )
  
  # calculate p values
  signif_data <- signif_data |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of(c('facet_id', sample_order))), groups) |>
    dplyr::summarise(
      test = list(test_fn(.data[[value]], ...)),
      p_value = test[[1]]$p.value,
      test = test[[1]]$method,
      .groups = 'drop'
    ) |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of(facet_cols))) |>
    dplyr::mutate(
      p_value = p.adjust(p_value, method = p_adjust),
      p_adjust_method = p_adjust
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::pick(tidyselect::any_of(c(facet_cols, sample_order))))
  
  # save test results
  if (! is.null(save_test_results)) {
    # select columns
    test_results <- signif_data |>
      dplyr::select(tidyselect::any_of(facet_cols), test, comparison = groups, p_value, tidyselect::any_of('p_adjust_method'))
    
    # save to `.pggf_test_results` environment
    assign(save_test_results, test_results |> dplyr::ungroup(), envir = .pggf_test_results)
  }
  
  # select final columns and split groups col
  signif_data <- signif_data |>
    tidyr::separate_wider_delim(
      groups,
      delim = ':',
      names = group
    ) |>
    dplyr::select(tidyselect::any_of(c(facet_cols, group)), p_value)
  
  return(signif_data)
}


# horizontal lines
pggf_hline <- function(pggf, fn, group = NULL, extra_cols = NULL) {
  # make sure `pggf` is a `pggf` class object
  .assert_pggf(pggf)
  
  # test if `fn` is given
  if (missing(fn) || .pggf_is_null(substitute(fn))) {
    cli::cli_abort('argument {.arg fn} is missing, with no default')
  }
  
  # calculate intercepts
  intercepts <- pggf$data |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'), !! pggf$facet_col, !! pggf$facet_row, !! pggf$facet_wrap, {{group}})) |>
    dplyr::summarise(
      intercept = {{fn}},
      dplyr::across(
        {{extra_cols}},
        ~ if (length(unique(.x)) == 1) {
          unique(.x)
        } else {
          cli::cli_warn('Non-unique values in {.arg extra_cols}')
          paste0(.x, collapse = ',')
        }
      ),
      .groups = 'drop'
    ) |>
    tidyr::unnest_longer(intercept)
  
  # export data
  intercept_file <- intercepts |>
    dplyr::mutate(
      file = paste0(pggf$filename, '_hline.tsv')
    ) |>
    dplyr::nest_by(file) |>
    dplyr::summarise(
      content = readr::format_tsv(data),
      .groups = 'drop'
    ) |>
    tibble::deframe()
  
  # add file to pggf object
  pggf$files <- c(pggf$files, intercept_file)
  
  # add to ggplot
  if (! is.null(pggf$plot)) {
    if (missing(group) || .pggf_is_null(substitute(group))) {
      pggf$plot <- append(pggf$plot, ggplot2::geom_hline(data = intercepts, mapping = ggplot2::aes(yintercept = intercept)))
    } else {
      plot_group <- .pggf_sym_or_str(substitute(group))
      pggf$plot <- append(pggf$plot, ggplot2::geom_hline(data = intercepts, mapping = ggplot2::aes(yintercept = intercept, color = {{plot_group}})))
    }
  }
  
  # add to code
  if (! is.null(pggf$code)) {
    if (missing(group) || .pggf_is_null(substitute(group))) {
      pggf$code <- append(pggf$code, '\t\\pggf_hline()\n', after = length(pggf$code) - 1)
    } else {
      group_str <- if (is.symbol(substitute(group))) {deparse(substitute(group))} else {as.character(group)}
      pggf$code <- append(pggf$code, paste0('\t\\pggf_hline(style from column = ', group_str, ')\n'), after = length(pggf$code) - 1) 
    }
  }
  
  # call the new `pggf` object
  pggf
}


# vertical lines
pggf_vline <- function(pggf, fn, group = NULL, extra_cols = NULL) {
  # make sure `pggf` is a `pggf` class object
  .assert_pggf(pggf)
  
  # test if `fn` is given
  if (missing(fn) || .pggf_is_null(substitute(fn))) {
    cli::cli_abort('argument {.arg fn} is missing, with no default')
  }
  
  # calculate intercepts
  intercepts <- pggf$data |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'), !! pggf$facet_col, !! pggf$facet_row, !! pggf$facet_wrap, {{group}})) |>
    dplyr::summarise(
      intercept = {{fn}},
      dplyr::across(
        {{extra_cols}},
        ~ if (length(unique(.x)) == 1) {
          unique(.x)
        } else {
          cli::cli_warn('Non-unique values in {.arg extra_cols}')
          paste0(.x, collapse = ',')
        }
      ),
      .groups = 'drop'
    ) |>
    tidyr::unnest_longer(intercept)
  
  # export data
  intercept_file <- intercepts |>
    dplyr::mutate(
      file = paste0(pggf$filename, '_vline.tsv')
    ) |>
    dplyr::nest_by(file) |>
    dplyr::summarise(
      content = readr::format_tsv(data),
      .groups = 'drop'
    ) |>
    tibble::deframe()
  
  # add file to pggf object
  pggf$files <- c(pggf$files, intercept_file)
  
  # add to ggplot
  if (! is.null(pggf$plot)) {
    if (missing(group) || .pggf_is_null(substitute(group))) {
      pggf$plot <- append(pggf$plot, ggplot2::geom_vline(data = intercepts, mapping = ggplot2::aes(xintercept = intercept)))
    } else {
      plot_group <- .pggf_sym_or_str(substitute(group))
      pggf$plot <- append(pggf$plot, ggplot2::geom_vline(data = intercepts, mapping = ggplot2::aes(xintercept = intercept, color = {{plot_group}})))
    }
  }
  
  # add to code
  if (! is.null(pggf$code)) {
    if (missing(group) || .pggf_is_null(substitute(group))) {
      pggf$code <- append(pggf$code, '\t\\pggf_vline()\n', after = length(pggf$code) - 1)
    } else {
      group_str <- if (is.symbol(substitute(group))) {deparse(substitute(group))} else {as.character(group)}
      pggf$code <- append(pggf$code, paste0('\t\\pggf_vline(style from column = ', group_str, ')\n'), after = length(pggf$code) - 1) 
    }
  }
  
  # call the new `pggf` object
  pggf
}


# linear regression line
pggf_trendline <- function(pggf, group = NULL, extra_cols = NULL, force_through = NULL) {
  # make sure `pggf` is a `pggf` class object
  .assert_pggf(pggf, dims = 2)
  
  # define model formula
  if (missing(force_through) || .pggf_is_null(substitute(force_through))) {
    lmform <- 'y ~ x'
    force_through <- NULL
  } else {
    if (! is.numeric(force_through) || length(force_through) != 1) {
      cli::cli_abort('{.arg force_through} must be a {.cls numeric} vector of length 1')
    }
    lmform <- paste0('I(y-', force_through, ') ~ 0+x')
  }
  
  # calculate slope and intercept
  trendline <- pggf$data |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'), !! pggf$facet_col, !! pggf$facet_row, !! pggf$facet_wrap, {{group}})) |>
    dplyr::summarise(
      slope = ifelse(is.null(force_through), coef(lm(as.formula(lmform)))[2], coef(lm(as.formula(lmform)))[1]),
      intercept = ifelse(is.null(force_through), coef(lm(as.formula(lmform)))[1], force_through),
      dplyr::across(
        {{extra_cols}},
        ~ if (length(unique(.x)) == 1) {
          unique(.x)
        } else {
          cli::cli_warn('Non-unique values in {.arg extra_cols}')
          paste0(.x, collapse = ',')
        }
      ),
      .groups = 'drop'
    )
  
  # export data
  trendline_file <- trendline |>
    dplyr::mutate(
      file = paste0(pggf$filename, '_trendline.tsv')
    ) |>
    dplyr::nest_by(file) |>
    dplyr::summarise(
      content = readr::format_tsv(data),
      .groups = 'drop'
    ) |>
    tibble::deframe()
  
  # add file to pggf object
  pggf$files <- c(pggf$files, trendline_file)
  
  # add to ggplot
  if (! is.null(pggf$plot)) {
    if (missing(group) || .pggf_is_null(substitute(group))) {
      pggf$plot <- append(pggf$plot, ggplot2::geom_smooth(method = 'lm', formula = lmform, se = FALSE))
    } else {
      plot_group <- .pggf_sym_or_str(substitute(group))
      pggf$plot <- append(pggf$plot, ggplot2::geom_smooth(mapping = ggplot2::aes(color = {{plot_group}}), method = 'lm', formula = lmform, se = FALSE))
    }
  }
  
  # add to code
  if (! is.null(pggf$code)) {
    if (missing(group) || .pggf_is_null(substitute(group))) {
      pggf$code <- append(pggf$code, '\t\\pggf_trendline()\n', after = length(pggf$code) - 1)
    } else {
      group_str <- if (is.symbol(substitute(group))) {deparse(substitute(group))} else {as.character(group)}
      pggf$code <- append(pggf$code, paste0('\t\\pggf_trendline(style from column = ', group_str, ')\n'), after = length(pggf$code) - 1) 
    }
  }
  
  # call the new `pggf` object
  pggf
}


# sloped lines
pggf_abline <- function(pggf, slope = 1, intercept = 0, group = NULL, extra_cols = NULL) {
  # make sure `pggf` is a `pggf` class object
  .assert_pggf(pggf)
  
  # calculate slope and/or intercept
  data_abline <- pggf$data |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'), !! pggf$facet_col, !! pggf$facet_row, !! pggf$facet_wrap, {{group}})) |>
    dplyr::summarise(
      slope = {{slope}},
      intercept = {{intercept}},
      dplyr::across(
        {{extra_cols}},
        ~ if (length(unique(.x)) == 1) {
          unique(.x)
        } else {
          cli::cli_warn('Non-unique values in {.arg extra_cols}')
          paste0(.x, collapse = ',')
        }
      ),
      .groups = 'drop'
    )
  
  # export data
  abline_file <- data_abline |>
    dplyr::mutate(
      file = paste0(pggf$filename, '_abline.tsv')
    ) |>
    dplyr::nest_by(file) |>
    dplyr::summarise(
      content = readr::format_tsv(data),
      .groups = 'drop'
    ) |>
    tibble::deframe()
  
  # add file to pggf object
  pggf$files <- c(pggf$files, abline_file)
  
  # add to ggplot
  if (! is.null(pggf$plot)) {
    if (missing(group) || .pggf_is_null(substitute(group))) {
      pggf$plot <- append(pggf$plot, ggplot2::geom_abline(data = data_abline, mapping = ggplot2::aes(slope = slope, intercept = intercept)))
    } else {
      plot_group <- .pggf_sym_or_str(substitute(group))
      pggf$plot <- append(pggf$plot, ggplot2::geom_abline(data = data_abline, mapping = ggplot2::aes(slope = slope, intercept = intercept, color = {{plot_group}})))
    }
  }
  
  # add to code
  if (! is.null(pggf$code)) {
    if (missing(group) || .pggf_is_null(substitute(group))) {
      pggf$code <- append(pggf$code, '\t\\pggf_abline()\n', after = length(pggf$code) - 1)
    } else {
      group_str <- if (is.symbol(substitute(group))) {deparse(substitute(group))} else {as.character(group)}
      pggf$code <- append(pggf$code, paste0('\t\\pggf_abline(style from column = ', group_str, ')\n'), after = length(pggf$code) - 1) 
    }
  }
  
  # call the new `pggf` object
  pggf
}


# quiver plot
pggf_quiver <- function(pggf, u = NULL, v = NULL, split = NULL, extra_cols = NULL) {
  # make sure `pggf` is a `pggf` class object
  .assert_pggf(pggf, dims = 2)
  
  # at least one of `u` and `v` is required
  if ((missing(u) || .pggf_is_null(substitute(u))) && (missing(v) || .pggf_is_null(substitute(v)))) {
    cli::cli_abort('At least one of {.arg u} or {.arg v} must be given')
  }
  
  # save `u`/`v` for axis limit calculations
  if (! missing(u) && ! .pggf_is_null(substitute(u))) {
    u_data <- pggf$data |>
      dplyr::select({{u}})
    
    if (ncol(u_data) == 1) {
      pggf$error[1] <- -u_data
      pggf$error[2] <- u_data
    } else {
      cli::cli_abort('{.arg u} must be of length 1')
    }
  }
  if (! missing(v) && ! .pggf_is_null(substitute(v))) {
    v_data <- pggf$data |>
      dplyr::select({{v}})
    
    if (ncol(v_data)  == 1) {
      pggf$error[3] <- -v_data
      pggf$error[4] <- v_data
    } else {
      cli::cli_abort('{.arg v} must be of length 1')
    }
  }
  
  # prepare `data`
  quiver_data <- pggf$data |>
    dplyr::select(tidyselect::any_of('facet_id'), x, y, 'u' = {{u}}, 'v' = {{v}}, {{split}}, {{extra_cols}})
  
  if (missing(u) || .pggf_is_null(substitute(u))) {
    quiver_data <- quiver_data |>
      dplyr::mutate(
        u = 0
      )
  }
  if (missing(v) || .pggf_is_null(substitute(v))) {
    quiver_data <- quiver_data |>
      dplyr::mutate(
        v = 0
      )
  }
  
  # split (optional) and export data
  files <- quiver_data |>
    dplyr::arrange(dplyr::pick(tidyselect::any_of('facet_id'))) |>
    dplyr::mutate('.file_basename' = pggf$filename, '.file_suffix' = 'quiver.tsv') |>
    tidyr::unite('file', .file_basename, {{split}}, .file_suffix, remove = FALSE) |>
    dplyr::select(-c(.file_basename, .file_suffix)) |>
    dplyr::nest_by(file) |>
    dplyr::summarise(
      content = readr::format_tsv(data),
      .groups = 'drop'
    ) |>
    tibble::deframe()
  
  # add file(s) to pggf object
  pggf$files <- c(pggf$files, files)
  
  # add to ggplot
  if (! is.null(pggf$plot)) {
    if (missing(u) || .pggf_is_null(substitute(u))) {
      plot_u <- 0
    } else {
      plot_u <- .pggf_sym_or_str(substitute(u))
    }
    if (missing(v) || .pggf_is_null(substitute(v))) {
      plot_v <- 0
    } else {
      plot_v <- .pggf_sym_or_str(substitute(v))
    }
    
    if (missing(split) || .pggf_is_null(substitute(split))) {
      pggf$plot <- append(
        pggf$plot,
        ggplot2::geom_segment(
          ggplot2::aes(xend = as.numeric(x) + {{plot_u}}, yend = as.numeric(y) + {{plot_v}}),
          arrow = ggplot2::arrow(length = ggplot2::unit(0.15, 'cm'))
        )
      )
    } else {
      plot_split <- .pggf_sym_or_str(substitute(split))
      pggf$plot <- append(
        pggf$plot,
        ggplot2::geom_segment(
          ggplot2::aes(xend = as.numeric(x) + {{plot_u}}, yend = as.numeric(y) + {{plot_v}}, color = {{plot_split}}),
          arrow = ggplot2::arrow(length = ggplot2::unit(0.15, 'cm'))
        )
      )
    }
  }
  
  # add to code
  if (! is.null(pggf$code)) {
    if (missing(split) || .pggf_is_null(substitute(split))) {
      pggf$code <- append(
        pggf$code,
        paste0('\t\\pggf_quiver()\n'),
        after = length(pggf$code) - 1
      )
    } else {
      code_data <- paste0('data = ', gsub('_quiver.tsv', '', names(files)))
      pggf$code <- append(
        pggf$code,
        c(paste0('\t\\pggf_quiver(', code_data, ')'), '\t'),
        after = length(pggf$code) - 1
      )
    }
  }
  
  # call the new `pggf` object
  pggf
}

# hexbin plot
pggf_hexbin <- function(pggf, bins = 50, count = 'log2', colorbar = 'fixed') {
  # make sure `pggf` is a `pggf` class object
  .assert_pggf(pggf, dims = 2)
  
  # helper for count functions
  count_fns <- c('raw' = identity, 'log2' = log2, 'log10' = log10)
  
  # check arguments
  if (! is.numeric(bins) || length(bins) != 1) {
    cli::cli_abort('{.arg bins} must be a {.cls numeric} vector of length 1')
  }
  if (! count %in% names(count_fns)) {
    cli::cli_abort('{.arg count} must be one of "raw", "log2", or "log10"')
  } else {
    count_fn <- count_fns[[count]]
  }
  if (! colorbar %in% c('fixed', 'free')) {
    cli::cli_abort('{.arg colorbar} must be "fixed" or "free"')
  }
  
  # get facet info
  if (! is.null(pggf$facets)) {
    has_cols <- ('col' %in% names(pggf$facets))
    has_rows <- ('row' %in% names(pggf$facets))
  } else {
    has_cols <- FALSE
    has_rows <- FALSE
  }
  
  # prepare `data`
  hexbin_data <- pggf$data |>
    dplyr::ungroup()
  
  # get x and y bounds
  if (pggf$scales == 'square') {
    hexbin_data <- hexbin_data |>
      dplyr::mutate(
        x_bounds = list(range(x, y)),
        y_bounds = x_bounds
      )
  } else {
    if (pggf$scales %in% c('free', 'free_x') && has_cols) {
      hexbin_data <- hexbin_data |>
        dplyr::group_by(!! pggf$facet_col) |>
        dplyr::mutate(
          x_bounds = list(range(x))
        ) |>
        dplyr::ungroup()
    } else {
      hexbin_data <- hexbin_data |>
        dplyr::mutate(
          x_bounds = list(range(x))
        )
    }
    if (pggf$scales %in% c('free', 'free_y') && has_rows) {
      hexbin_data <- hexbin_data |>
        dplyr::group_by(!! pggf$facet_row) |>
        dplyr::mutate(
          y_bounds = list(range(y))
        ) |>
        dplyr::ungroup()
    } else {
      hexbin_data <- hexbin_data |>
        dplyr::mutate(
          y_bounds = list(range(y))
        )
    }
  }
  
  # get hexbin coordinates
  hexbin_data <- hexbin_data |>
    dplyr::nest_by(dplyr::pick(tidyselect::any_of('facet_id'), !! pggf$facet_col, !! pggf$facet_row)) |>
    dplyr::reframe(
      hexbin = list(
        hexbin::hexbin(
          x = data$x,
          y = data$y,
          xbins = bins,
          xbnds = range(data$x_bounds),
          ybnds = range(data$y_bounds)
        )
      ),
      tibble::as_tibble(hexbin::hcell2xy(hexbin)),
      count = slot(hexbin, 'count'),
      count = count_fn(count),
      .pggf.dim.x = slot(hexbin, 'dimen')[2],
      .pggf.dim.y = slot(hexbin, 'dimen')[1],
      .pggf.hex.width = (max(slot(hexbin, 'xbnds')) - min(slot(hexbin, 'xbnds'))) / (.pggf.dim.x - 1),
      .pggf.hex.shape = diff(slot(hexbin, 'ybnds')) / diff(slot(hexbin, 'xbnds'))
    ) |>
    dplyr::select(-hexbin)
  
  # export data
  hexbin_file <- hexbin_data |>
    dplyr::select(-tidyselect::starts_with('.pggf.')) |>
    dplyr::mutate(
      file = paste0(pggf$filename, '_hexbin.tsv')
    ) |>
    dplyr::nest_by(file) |>
    dplyr::summarise(
      content = readr::format_tsv(data),
      .groups = 'drop'
    ) |>
    tibble::deframe()
  
  # add file to pggf object
  pggf$files <- c(pggf$files, hexbin_file)
  
  # calculate colorbar limits
  colorbar_limits <- hexbin_data |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
    dplyr::summarise(
      dplyr::across(count, list('min' = ~ min(.x, na.rm = TRUE), 'max' = ~ max(.x, na.rm = TRUE)), .names = 'point meta {.fn}'),
      dplyr::across(tidyselect::starts_with('.pggf.'), unique),
      .groups = 'drop'
    ) |>
    dplyr::ungroup()
  
  # add all possible facets
  if (! is.null(pggf$facets) && ! 'wrap' %in% names(pggf$facets)) {
    colorbar_limits <- colorbar_limits |>
      dplyr::full_join(dplyr::select(pggf$facets, facet_id), dplyr::join_by(facet_id))
  }
  
  # set global colorbar limits for 'fixed' colorbars
  if (colorbar == 'fixed') {
    colorbar_limits <- colorbar_limits |>
      dplyr::mutate(
        `point meta min` = min(`point meta min`, na.rm = TRUE),
        `point meta max` = max(`point meta max`, na.rm = TRUE)
      )
  } else { # make sure there are no NAs in the colorbar limits (e.g. from empty facets)
    colorbar_limits <- colorbar_limits |>
      dplyr::mutate(
        `point meta min` = tidyr::replace_na(`point meta min`, 0),
        `point meta max` = tidyr::replace_na(`point meta max`, 1)
      )
  }
  
  # export colorbar limits  
  pggf$axes_extra <- colorbar_limits
  
  # add to ggplot
  if (! is.null(pggf$plot)) {
    pggf$plot <- append(
      pggf$plot,
      c(
        ggplot2::geom_hex(bins = bins),
        ggplot2::scale_fill_viridis_c(trans = ifelse(count == 'raw', 'identity', count))
      )
    )
  }
  
  # add to code
  if (! is.null(pggf$code)) {
    pggf$code <- append(
      pggf$code,
      paste0('\t', ifelse(colorbar == 'fixed', 'facet 1', 'every facet'), '/.append style = {pggf colorbar = {title = count, ', count,'}},'),
      after = 1
    )
    pggf$code <- append(
      pggf$code,
      paste0('\t\\pggf_hexbin(bins = ', bins, ')\n'),
      after = length(pggf$code) - 1
    )
  }
  
  # call the new `pggf` object
  pggf
}

# bar plot
pggf_bar <- function(pggf, group = NULL, error = NULL, extra_cols = NULL, orientation = 'x', stack = FALSE) {
  # make sure `pggf` is a `pggf` class object
  .assert_pggf(pggf, dims = 2)
  
  # throw an error if `error` is used together with `stack = TRUE`
  if (stack && (! missing(error) && ! .pggf_is_null(substitute(error)))) {
    cli::cli_abort('Error bars cannot be added to stacked bars')
  }
  
  # get orientation
  if (orientation %in% c('x', 'v', 'vertical')) {
    orientation <- 'x'
    if (! pggf$numeric_y) {
      cli::cli_abort('{.arg y} must be numeric for {.arg orientation = "x"}')
    }
    sample_cols <- c('x', 'xticklabels')
    value_col <- 'y'
  } else if (orientation %in% c('y', 'h', 'horizontal')) {
    orientation <- 'y'
    if (! pggf$numeric_x) {
      cli::cli_abort('{.arg x} must be numeric for {.arg orientation = "y"}')
    }
    sample_cols <- c('y', 'yticklabels')
    value_col <- 'x'
  } else {
    cli::cli_abort(c('{.arg orientation} must be "x" or "y"', i = 'You can also use "v" or "vertical" for "x" and "h" or "horizontal" for "y"'))
  }
  
  # save error data for axis limit calculations
  if (! missing(error) && ! .pggf_is_null(substitute(error))) {
    error_data <- pggf$data |>
      dplyr::select({{error}})
    
    if (ncol(error_data) %in% c(1, 2)) {
      if (orientation == 'x') {
        pggf$error[c(3, 4)] <- error_data
      } else {
        pggf$error[c(1, 2)] <- error_data
      }
    } else {
      cli::cli_abort('{.arg error} must be of length 1 or 2')
    }
  }
  
  # save stack height for axis limit calculations
  if (stack) {
    pggf$extra_coords <- pggf$extra_coords |>
      dplyr::bind_rows(
        pggf$data |>
          dplyr::group_by(
            dplyr::pick(
              tidyselect::any_of(c('facet_id', orientation, paste0(orientation, 'ticklabels'))),
              !! pggf$facet_col,
              !! pggf$facet_row,
              !! pggf$facet_wrap
            )
          ) |>
          dplyr::summarise(
            dplyr::across(tidyselect::all_of(value_col), sum),
            .groups = 'drop'
          )
      )
  }
  
  # prepare `data`
  bar_data <- pggf$data |>
    dplyr::select(tidyselect::any_of(c('facet_id', sample_cols, value_col)), 'group' = {{group}}, {{error}}, {{extra_cols}}) |>
    dplyr::mutate(
      dplyr::across(tidyselect::any_of('group'), as.ordered)
    ) |>
    dplyr::arrange(dplyr::pick(tidyselect::any_of(c('facet_id', 'group', orientation))))
  
  # calculate bar width and shift
  shift_col <- paste0(orientation, '_shift')
  
  if (stack) {
    bar_data <- bar_data |>
      dplyr::mutate(
        bar_width = 1,
        {{shift_col}} := 0
      )
  } else {
    bar_data <- bar_data |>
      dplyr::mutate(
        group_n = 1,
        {{shift_col}} := 0,
        dplyr::across(tidyselect::any_of('group'), ~ max(as.numeric(droplevels(.x))), .names = 'group_n'),
        bar_width = 1 / group_n,
        dplyr::across(tidyselect::any_of('group'), ~ (as.numeric(droplevels(.x)) - 0.5 * (group_n + 1)) * bar_width, .names = shift_col),
      ) |>
      dplyr::select(-group_n)
  }
  
  # export data
  bar_file <- bar_data |>
    dplyr::mutate(
      dplyr::across(tidyselect::ends_with('ticklabels'), ~ dplyr::if_else(grepl(' ', .x, fixed = TRUE), paste0('{', .x, '}'), .x)),
      file = paste0(pggf$filename, '_bar.tsv')
    ) |>
    dplyr::nest_by(file) |>
    dplyr::summarise(
      content = readr::format_tsv(data),
      .groups = 'drop'
    ) |>
    tibble::deframe()
  
  # add file to pggf object
  pggf$files <- c(pggf$files, bar_file)
  
  # add to ggplot
  if (! is.null(pggf$plot)) {
    if (! missing(group) && ! .pggf_is_null(substitute(group))) {
      plot_group <- .pggf_sym_or_str(substitute(group))
      pggf$plot <- append(
        pggf$plot,
        ggplot2::geom_col(
          mapping = ggplot2::aes(fill = as.ordered({{plot_group}})),
          position = ifelse(stack, 'stack', 'dodge'),
          orientation = orientation
        )
      )
    } else {
      pggf$plot <- append(
        pggf$plot,
        ggplot2::geom_col(
          position = ifelse(stack, 'stack', 'dodge'),
          orientation = orientation
        )
      )
    }
    if (! missing(error) && ! .pggf_is_null(substitute(error))) {
      plot_err <- names(error_data)
      if (length(plot_err) == 1) {
        plot_err[2] <- plot_err
      }
      if (orientation == 'x') {
        if (! missing(group) && ! .pggf_is_null(substitute(group))) {
          pggf$plot <- append(
            pggf$plot,
            ggplot2::geom_errorbar(
              ggplot2::aes(ymin = y - .data[[plot_err[1]]], ymax = y + .data[[plot_err[2]]], color = as.ordered({{plot_group}})),
              position = ggplot2::position_dodge(0.9),
              width = 0.1
            )
          )
        } else {
          pggf$plot <- append(
            pggf$plot,
            ggplot2::geom_errorbar(
              ggplot2::aes(ymin = y - .data[[plot_err[1]]], ymax = y + .data[[plot_err[2]]]),
              width = 0.1
            )
          )
        }
      } else {
        if (! missing(group) && ! .pggf_is_null(substitute(group))) {
          pggf$plot <- append(
            pggf$plot,
            ggplot2::geom_errorbar(
              ggplot2::aes(xmin = x - .data[[plot_err[1]]], xmax = x + .data[[plot_err[2]]], color = as.ordered({{plot_group}})),
              position = ggplot2::position_dodge(0.9),
              width = 0.1
            )
          )
        } else {
          pggf$plot <- append(
            pggf$plot,
            ggplot2::geom_errorbar(
              ggplot2::aes(xmin = x - .data[[plot_err[1]]], xmax = x + .data[[plot_err[2]]]),
              width = 0.1
            )
          )
        }
      }
    }
  }
  
  # add to code
  if (! is.null(pggf$code)) {
    if (! missing(error) && ! .pggf_is_null(substitute(error))) {
      if (length(names(error_data)) == 1) {
        code_err <- paste0(', ', value_col, ' error = {', names(error_data), '}')
      } else {
        code_err <- paste0(', ', value_col, ' error* = {', names(error_data)[1], '}{', names(error_data)[2], '}')
      }
    } else {
      code_err <- NULL
    }
    
    pggf$code <- append(
      pggf$code,
      paste0('\t\\pggf_bar(', value_col, 'bar', ifelse(stack, ' stacked', ''), code_err, ')\n'),
      after = length(pggf$code) - 1
    )
  }
  
  # call the new `pggf` object
  pggf
}


# heatmap
pggf_heatmap <- function(pggf, score, extra_cols = NULL, colorbar = 'fixed', dendrogram = 'none') {
  # make sure `pggf` is a `pggf` class object
  .assert_pggf(pggf, dims = 2)
  
  # make sure `x` and `y` are non-numeric
  if (pggf$numeric_x || pggf$numeric_y) {
    cli::cli_abort('{.fn pggf_heatmap} only works with non-numeric x and y data')
  }
  
  # check arguments
  if (missing(score)) {
    cli::cli_abort('Required argument {.arg score} is missing')
  }
  if (! colorbar %in% c('fixed', 'free')) {
    cli::cli_abort('{.arg colorbar} must be "fixed" or "free"')
  }
  else if (! dendrogram %in% c('none', 'x', 'y', 'both')) {
    cli::cli_abort('{.arg dendrogram} must be one of: "none", "x", "y", or "both"')
  }
  
  # get facet info
  if (! is.null(pggf$facets)) {
    has_cols <- ('col' %in% names(pggf$facets))
    has_rows <- ('row' %in% names(pggf$facets))
  } else {
    has_cols <- FALSE
    has_rows <- FALSE
  }
  
  # prepare `data`
  heatmap_data <- pggf$data |>
    dplyr::ungroup()
  
  # add all possible x/y combinations
  if (is.null(pggf$facets)) {
    all_combis <- heatmap_data |>
      tidyr::expand(x, y)
    heatmap_data <- all_combis |>
      dplyr::left_join(
        heatmap_data |>
          dplyr::select(x, xticklabels, y, yticklabels, {{score}}, {{extra_cols}}),
        dplyr::join_by(x, y)
      )
  } else {
    if (has_cols && pggf$scales %in% c('free', 'free_x')) {
      all_combis <- heatmap_data |>
        dplyr::distinct(!! pggf$facet_col, x) |>
        dplyr::full_join(
          pggf$facets,
          dplyr::join_by(!! pggf$facet_col == col),
          relationship = 'many-to-many'
        )
    } else {
      all_combis <- heatmap_data |>
        dplyr::distinct(x) |>
        tidyr::expand_grid(pggf$facets)
    }
    if (has_rows && pggf$scales %in% c('free', 'free_y')) {
      all_combis <- heatmap_data |>
        dplyr::distinct(!! pggf$facet_row, y) |>
        dplyr::full_join(
          all_combis,
          dplyr::join_by(!! pggf$facet_row == row),
          relationship = 'many-to-many'
        )
    } else {
      all_combis <- heatmap_data |>
        dplyr::distinct(y) |>
        tidyr::expand_grid(all_combis)
    }
    heatmap_data <- all_combis |>
      dplyr::select(facet_id, x, y) |>
      dplyr::left_join(
        heatmap_data |>
          dplyr::select(facet_id, x, xticklabels, y, yticklabels, {{score}}, {{extra_cols}}),
        dplyr::join_by(facet_id, x, y)
      )
  }
  
  # create dendrograms
  x_dendrogram <- NULL
  y_dendrogram <- NULL
  
  if (dendrogram %in% c('x', 'both')) {
    if (is.null(pggf$facets) || (! has_rows && pggf$scales %in% c('free', 'free_x'))) {
      heatmap_data <- heatmap_data |>
        dplyr::group_by(!! pggf$facet_row) |>
        dplyr::mutate(
          .pggf_dendrogram(xticklabels, y, {{score}}),
        ) |>
        dplyr::mutate(
          x = as.numeric(ordered(xticklabels, levels = order[[1]]))
        ) |>
        dplyr::ungroup()
      
      x_dendrogram <- heatmap_data |>
        dplyr::distinct(!! pggf$facet_row, 'x dendrogram' = `axis option`)
      
      heatmap_data <- heatmap_data |>
        dplyr::select(-c(order, `axis option`))
      
      # replace original data with reordered data
      pggf$data <- heatmap_data
    } else {
      cli::cli_abort('Dendrograms for {.arg x} can only be used when there are no facet columns and {.arg scales = "free"} or {.arg scales = "free_x"} was used')
    }
  }
  
  if (dendrogram %in% c('y', 'both')) {
    if (is.null(pggf$facets) || (! has_rows && pggf$scales %in% c('free', 'free_x'))) {
      heatmap_data <- heatmap_data |>
        dplyr::group_by(!! pggf$facet_col) |>
        dplyr::mutate(
          .pggf_dendrogram(yticklabels, x, {{score}}),
        ) |>
        dplyr::mutate(
          y = as.numeric(ordered(yticklabels, levels = order[[1]]))
        ) |>
        dplyr::ungroup()
      
      y_dendrogram <- heatmap_data |>
        dplyr::distinct(dplyr::pick(tidyselect::any_of('facet_id')), 'y dendrogram' = `axis option`)
      
      heatmap_data <- heatmap_data |>
        dplyr::select(-c(order, `axis option`))
      
      # replace original data with reordered data
      pggf$data <- heatmap_data
    } else {
      cli::cli_abort('Dendrograms for {.arg y} can only be used when there are no facet rows and {.arg scales = "free"} or {.arg scales = "free_y"} was used')
    }
  }
  
  # remove ticklabels and order data
  heatmap_data <- heatmap_data |>
    dplyr::select(-xticklabels, -yticklabels) |>
    dplyr::rename('score' = {{score}}) |>
    dplyr::arrange(dplyr::pick(tidyselect::any_of('facet_id')), x, y)
  
  # export data
  heatmap_file <- heatmap_data |>
    dplyr::mutate(
      file = paste0(pggf$filename, '_heatmap.tsv')
    ) |>
    dplyr::nest_by(file) |>
    dplyr::summarise(
      content = readr::format_tsv(data, na = 'NaN'),
      .groups = 'drop'
    ) |>
    tibble::deframe()
  
  # add file to pggf object
  pggf$files <- c(pggf$files, heatmap_file)
  
  # export column and row numbers, colorbar limits, and dendrograms
  heatmap_axis <- heatmap_data |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
    dplyr::summarise(
      `mesh/cols` = dplyr::n_distinct(x),
      `mesh/rows` = dplyr::n_distinct(y),
      dplyr::across(score, list('min' = ~ min(.x, na.rm = TRUE), 'max' = ~ max(.x, na.rm = TRUE)), .names = 'point meta {.fn}'),
      .groups = 'drop'
    ) |>
    dplyr::ungroup()
  
  if (colorbar == 'fixed') {
    heatmap_axis <- heatmap_axis|>
      dplyr::mutate(
        `point meta min` = min(`point meta min`, na.rm = TRUE),
        `point meta max` = max(`point meta max`, na.rm = TRUE)
      )
  }
  
  if (! is.null(x_dendrogram)) {
    if ('facet_id' %in% names(heatmap_axis)) {
      heatmap_axis <- heatmap_axis |>
        dplyr::inner_join(
          x_dendrogram,
          by = 'facet_id'
        )
    } else {
     heatmap_axis <- heatmap_axis |>
       dplyr::bind_cols(x_dendrogram)
    }
  }
  if (! is.null(y_dendrogram)) {
    if ('facet_id' %in% names(heatmap_axis)) {
      heatmap_axis <- heatmap_axis |>
        dplyr::inner_join(
          y_dendrogram,
          by = 'facet_id'
        )
    } else {
      heatmap_axis <- heatmap_axis |>
        dplyr::bind_cols(y_dendrogram)
    }
  }
  
  pggf$axes_extra <- heatmap_axis |>
    dplyr::mutate(
      `colormap limits` = paste(`point meta min`, `point meta max`, sep = ',')
    )
  
  # add to ggplot
  if (! is.null(pggf$plot)) {
    plot_fill <- .pggf_sym_or_str(substitute(score))
    pggf$plot <- append(pggf$plot, ggplot2::geom_tile(ggplot2::aes(fill = {{plot_fill}})))
    if (pggf$numeric_y) {
      pggf$plot <- append(pggf$plot, ggplot2::scale_y_continuous(limits = rev))
    } else {
      pggf$plot <- append(pggf$plot, ggplot2::scale_y_discrete(limits = rev))
    }
  }
  
  # add to code
  if (! is.null(pggf$code)) {
    code_dendrogram <- NULL
    if (! is.null(x_dendrogram)) {
      code_dendrogram <- c(code_dendrogram, '\tshow x dendrogram,')
    }
    if (! is.null(y_dendrogram)) {
      code_dendrogram <- c(code_dendrogram, '\tshow y dendrogram,')
    }
    
    code_title <- if (is.symbol(substitute(score))) {deparse(substitute(score))} else {as.character(score)}
    
    pggf$code <- append(
      pggf$code,
      c(
        code_dendrogram,
        paste0('\t', ifelse(colorbar == 'fixed', 'facet 1', 'every facet'), '/.append style = {pggf colorbar = {title = ', code_title, '}},')
      ),
      after = 1
    )
    pggf$code <- append(
      pggf$code,
      paste0('\t\\pggf_heatmap()\n'),
      after = length(pggf$code) - 1
    )
  }
  
  # call the new `pggf` object
  pggf
}


# helper function to create dendrograms
.pggf_dendrogram <- function(x, y, z) {
  dendrogram <- tibble::tibble(x = x, y = y, z = z) |>
    tidyr::pivot_wider(
      names_from = x,
      values_from = z
    ) |>
    dplyr::select(-y) |>
    as.matrix() |>
    t() |>
    dist() |>
    hclust()
  
  dendrogram_coords <- dendrogram$merge |>
    tibble::as_tibble(.name_repair = ~ paste0('cluster', seq_len(2))) |>
    dplyr::mutate(
      across(starts_with('cluster'), ~ dplyr::if_else(.x < 0, paste0('l', sapply(abs(.x), function(x) which(x == dendrogram$order) - 1)), paste0('c', .x))),
      height = dendrogram$height / max(dendrogram$height)
    ) |>
    dplyr::summarise(
      coords = paste(cluster1, cluster2, height, sep = '/', collapse = ','),
      .groups = 'drop'
    ) |>
    dplyr::pull(coords)
  
  dendrogram_order <- dendrogram$labels[dendrogram$order]
  
  return(tibble::tibble('axis option' = dendrogram_coords, 'order' = list(dendrogram_order)))
}


# line plot
pggf_line <- function(pggf, group = NULL, split = NULL, extra_cols = NULL, error = NULL, orientation = 'x') {
  # make sure `pggf` is a `pggf` class object
  .assert_pggf(pggf, dims = 2)
  
  # get orientation
  if (orientation %in% c('x', 'v', 'vertical')) {
    orientation <- 'x'
    if (! pggf$numeric_y) {
      cli::cli_abort('{.arg y} must be numeric for {.arg orientation = "x"}')
    }
    value_col <- 'y'
  } else if (orientation %in% c('y', 'h', 'horizontal')) {
    orientation <- 'y'
    if (! pggf$numeric_x) {
      cli::cli_abort('{.arg x} must be numeric for {.arg orientation = "y"}')
    }
    value_col <- 'x'
  } else {
    cli::cli_abort(c('{.arg orientation} must be "x" or "y"', i = 'You can also use "v" or "vertical" for "x" and "h" or "horizontal" for "y"'))
  }
  
  # prepare `data`
  data <- pggf$data |>
    dplyr::ungroup() |>
    dplyr::select(tidyselect::any_of('facet_id'), x, y, {{group}}, {{split}}, {{extra_cols}}, {{error}})
  
  # save error data for axis limit calculations and calculate error coordinates
  if (! missing(error) && ! .pggf_is_null(substitute(error))) {
    error_data <- data |>
      dplyr::select({{error}})
    
    if (ncol(error_data) %in% c(1, 2)) {
      if (orientation == 'x') {
        pggf$error[c(3, 4)] <- error_data
      } else {
        pggf$error[c(1, 2)] <- error_data
      }
    } else {
      cli::cli_abort('{.arg error} must be of length 1 or 2')
    }
    
    error_names <- names(error_data)
    if (length(error_names) == 1) {
      error_names[2] <- error_names
    }
    
    data <- data |>
      dplyr::mutate(
        dplyr::across(
          tidyselect::all_of(value_col),
          list('min' = ~ .x - .data[[error_names[1]]], 'max' = ~ .x + .data[[error_names[2]]]),
          .names = 'error_{value_col}_{fn}'
        )
      )
    
    if (length(error_names) == 1) {
      data <- data |>
        dplyr::rename(
          !! paste0('error_', value_col) := tidyselect::all_of(error_names[1])
        )
    } else {
      data <- data |>
        dplyr::rename(
          !! paste0('error_', value_col, '_1') := tidyselect::all_of(error_names[1]),
          !! paste0('error_', value_col, '_2') := tidyselect::all_of(error_names[2]),
        )
    }
  }
  
  # order data
  data <- data |>
    dplyr::arrange(dplyr::pick(tidyselect::any_of('facet_id'), {{group}}, {{split}}, tidyselect::any_of(orientation)))
  
  # split (optional) and export data
  files <- data |>
    dplyr::mutate('.file_basename' = pggf$filename, '.file_suffix' = 'line.tsv') |>
    tidyr::unite('file', .file_basename, {{split}}, .file_suffix, remove = FALSE) |>
    dplyr::select(-c(.file_basename, .file_suffix)) |>
    dplyr::nest_by(file) |>
    dplyr::summarise(
      content = readr::format_tsv(data),
      .groups = 'drop'
    ) |>
    tibble::deframe()
  
  # add file(s) to pggf object
  pggf$files <- c(pggf$files, files)
  
  # add to ggplot
  if (! is.null(pggf$plot)) {
    plot_group <- .pggf_sym_or_str(substitute(group))
    pggf$plot <- append(pggf$plot, ggplot2::geom_line(mapping = ggplot2::aes(color = {{plot_group}}, group = {{plot_group}}), orientation = orientation))
    if (! missing(error) && ! .pggf_is_null(substitute(error))) {
      if (orientation == 'x') {
        pggf$plot <- append(pggf$plot, ggplot2::geom_line(ggplot2::aes(y = y - .data[[error_names[1]]], color = {{plot_group}}, group = {{plot_group}}), orientation = orientation, linetype = 'dashed'))
        pggf$plot <- append(pggf$plot, ggplot2::geom_line(ggplot2::aes(y = y + .data[[error_names[2]]], color = {{plot_group}}, group = {{plot_group}}), orientation = orientation, linetype = 'dashed'))
      } else {
        pggf$plot <- append(pggf$plot, ggplot2::geom_line(ggplot2::aes(x = y - .data[[error_names[1]]], color = {{plot_group}}, group = {{plot_group}}), orientation = orientation, linetype = 'dashed'))
        pggf$plot <- append(pggf$plot, ggplot2::geom_line(ggplot2::aes(x = y + .data[[error_names[2]]], color = {{plot_group}}, group = {{plot_group}}), orientation = orientation, linetype = 'dashed'))
      }
    }
  }
  
  # add to code
  if (! is.null(pggf$code)) {
    if (! missing(error) && ! .pggf_is_null(substitute(error))) {
      code_error <- paste0('error = lines')
    } else {
      code_error <- NULL
    }
    
    if (! missing(group) && ! .pggf_is_null(substitute(group))) {
      groups <- data |> dplyr::distinct(dplyr::pick({{group}})) |> dplyr::pull() |> sort() |> paste0(collapse = ',')
      group_str <- if (is.symbol(substitute(group))) {deparse(substitute(group))} else {as.character(group)}
      code_group <- paste0('group column = ', group_str)
      pggf$code <- append(
        pggf$code,
        paste0('\tline styles = {', groups, '},'),
        after = 1
      )
    } else {
      code_group <- NULL
    }
    
    code_extras <- paste0(c(code_error, code_group), collapse = ', ')
    
    if (! missing(split) && ! .pggf_is_null(substitute(split))) {
      if (code_extras != '') {
        code_extras <- paste0(', ', code_extras)
      }
      code_split <- paste0('data = ', gsub('_line.tsv', '', names(files)), code_extras)
      pggf$code <- append(
        pggf$code,
        c(paste0('\t\\pggf_line(', code_split, ')'), '\t'),
        after = length(pggf$code) - 1
      )
    } else {
      pggf$code <- append(
        pggf$code,
        paste0('\t\\pggf_line(', code_extras, ')\n'),
        after = length(pggf$code) - 1
      )
    }
    
    if (! is.null(code_group)) {
      pggf$code <- append(pggf$code, paste0('\t\\pggf_annotate(line legend)\n'), after = length(pggf$code) - 1)
    }
  }
  
  # call the new `pggf` object
  pggf
}


# density plots
pggf_density <- function(pggf, group = NULL, subdensity = pggf_defaults$subdensity,
                         stat = 'density', stack = FALSE, extra_cols = NULL,
                         samples = pggf_defaults$samples) {
  # make sure `pggf` is a `pggf` class object
  .assert_pggf(pggf, dims = 1)
  
  # check arguments
  if (! stat %in% c('density', 'count')) {
    cli::cli_abort('{.arg stat} must be "density" or "count"')
  }
  if (! is.numeric(samples) || length(samples) != 1) {
    cli::cli_abort('{.arg samples} must be a single number.')
  }
  
  # get orientation
  if (is.null(pggf$x_str)) {
    orientation = 'x'
    value_col = 'y'
  } else {
    orientation = 'y'
    value_col = 'x'
  }
  
  # make sure input data is numeric
  if (! pggf[[paste0('numeric_', value_col)]]) {
    cli::cli_abort('Numeric data is requiered for {.fn pggf_density}')
  }
  
  # calculate density
  density_data <- pggf$data |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
    dplyr::mutate(
      .pggf.bw = bw.nrd0(.data[[value_col]]),
      dplyr::across(
        tidyselect::all_of(value_col),
        list('min' = ~ min(.x - 3 * .pggf.bw), 'max' = ~ max(.x + 3 * .pggf.bw), 'n' = length),
        .names = '.pggf.{.fn}'
      )
    ) |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'), !! pggf$facet_col, !! pggf$facet_row, !! pggf$facet_wrap, {{group}})) |>
    dplyr::mutate(
      dplyr::across(
        {{extra_cols}},
        ~ if (length(unique(.x)) == 1) {
          unique(.x)
        } else {
          cli::cli_warn('Non-unique values in {.arg extra_cols}')
          paste0(.x, collapse = ',')
        }
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::nest_by(dplyr::pick(tidyselect::any_of('facet_id'), !! pggf$facet_col, !! pggf$facet_row, !! pggf$facet_wrap, {{group}}, {{extra_cols}})) |>
    dplyr::reframe(
      density = list(density(data[[value_col]], from = data$.pggf.min[1], to = data$.pggf.max[1], n = {{samples}})),
      x = density[[value_col]],
      y = density[[orientation]],
      n = length(data[[value_col]]),
      total_n = data$.pggf.n[1]
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-density)
  
  # scale densities
  if (stat == 'count') {
    density_data <- density_data |>
      dplyr::mutate(
        dplyr::across(
          tidyselect::all_of(orientation),
          ~ .x * n
        )
      )
  } else if (subdensity) {
    # normally, the density of each group has an area of 1
    # with `subdensity = TRUE`, the sum of the area of all group densities is 1
    # i.e. the subdensity takes into account the relative proportion of a group in the total data
    density_data <- density_data |>
      dplyr::mutate(
        dplyr::across(
          tidyselect::all_of(orientation),
          ~ .x * n / total_n
        )
      )
  }
  
  # select relevant columns and order data
  density_data <- density_data |>
    dplyr::select(tidyselect::any_of('facet_id'), !! pggf$facet_col, !! pggf$facet_row, !! pggf$facet_wrap, {{group}}, {{extra_cols}}, x, y) |>
    dplyr::arrange(dplyr::pick(tidyselect::any_of('facet_id')), {{group}})
  
  # export data
  density_file <- density_data |>
    dplyr::select(tidyselect::any_of('facet_id'), {{group}}, {{extra_cols}}, x, y) |>
    dplyr::mutate(
      file = paste0(pggf$filename, '_density.tsv')
    ) |>
    dplyr::nest_by(file) |>
    dplyr::summarise(
      content = readr::format_tsv(data),
      .groups = 'drop'
    ) |>
    tibble::deframe()
  
  # add file to pggf object
  pggf$files <- c(pggf$files, density_file)
  
  # replace original data with density data
  pggf$data <- density_data
  pggf[paste0(orientation, '_str')] <- stat
  
  # save stack height for axis limit calculations
  if (stack) {
    pggf$extra_coords <- pggf$extra_coords |>
      dplyr::bind_rows(
        pggf$data |>
          dplyr::group_by(
            dplyr::pick(
              tidyselect::any_of(c('facet_id', value_col, paste0(value_col, 'ticklabels'))),
              !! pggf$facet_col,
              !! pggf$facet_row,
              !! pggf$facet_wrap
            )
          ) |>
          dplyr::summarise(
            dplyr::across(tidyselect::all_of(orientation), sum),
            .groups = 'drop'
          )
      )
  }
  
  # add to ggplot
  if (! is.null(pggf$plot)) {
    plot_group <- .pggf_sym_or_str(substitute(group))
    pggf$plot <- append(
      pggf$plot,
      ggplot2::geom_area(
        mapping = ggplot2::aes(fill = {{plot_group}}, group = {{plot_group}}),
        orientation = value_col,
        position = dplyr::if_else(stack, 'stack', 'identity'),
        alpha = 0.5,
        color = 'black'
      )
    )
  }
  
  # add to code
  if (! is.null(pggf$code)) {
    if (! missing(group) && ! .pggf_is_null(substitute(group))) {
      groups <- density_data |> dplyr::distinct(dplyr::pick({{group}})) |> dplyr::pull() |> sort() |> paste0(collapse = ',')
      group_str <- if (is.symbol(substitute(group))) {deparse(substitute(group))} else {as.character(group)}
      code_group <- paste0('group column = ', group_str)
      pggf$code <- append(
        pggf$code,
        paste0('\tdensity styles = {', groups, '},'),
        after = 1
      )
    } else {
      code_group <- NULL
    }
    
    if (stack) {
      code_stack <- paste0('stack plots = ', orientation)
    } else {
      code_stack <- NULL
    }
    
    code_extras <- paste0(c(code_group, code_stack), collapse = ', ')
    
    pggf$code <- append(
      pggf$code,
      paste0('\t\\pggf_density(', code_extras, ')\n'),
      after = length(pggf$code) - 1
    )
    
    if (! is.null(code_group)) {
      pggf$code <- append(pggf$code, paste0('\t\\pggf_annotate(density legend)\n'), after = length(pggf$code) - 1)
    }
  }
  
  # call the new `pggf` object
  pggf
}


# ridgeline plots
pggf_ridgeline <- function(pggf, orientation = 'x', subdensity = pggf_defaults$subdensity, 
                           tails = 'extend', extra_cols = NULL, samples = pggf_defaults$samples) {
  # make sure `pggf` is a `pggf` class object
  .assert_pggf(pggf, dims = 2)
  
  # check arguments
  if (! is.numeric(samples) || length(samples) != 1) {
    cli::cli_abort('{.arg samples} must be a single number.')
  }
  if (! tails %in% c('extend', 'sample')) {
    cli::cli_abort('{.arg tails} must be one of "extend" or "sample"')
  }
  
  # get orientation
  if (orientation %in% c('x', 'h', 'horizontal')) {
    orientation <- 'x'
    if (! pggf$numeric_x) {
      cli::cli_abort('{.arg x} must be numeric for {.arg orientation = "x"}')
    }
    sample_cols <- c('y', 'yticklabels')
    value_col <- 'x'
  } else if (orientation %in% c('y', 'v', 'vertical')) {
    orientation <- 'y'
    if (! pggf$numeric_y) {
      cli::cli_abort('{.arg y} must be numeric for {.arg orientation = "y"}')
    }
    sample_cols <- c('x', 'xticklabels')
    value_col <- 'y'
  } else {
    cli::cli_abort(c('{.arg orientation} must be "x" or "y"', i = 'You can also use "h" or "horizontal" for "x" and "v" or "vertical" for "y"'))
  }
  
  # calculate densities
  ridgeline_data <- pggf$data |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
    dplyr::mutate(
      .pggf.bw = bw.nrd0(.data[[value_col]]),
      dplyr::across(
        tidyselect::all_of(value_col),
        list('min' = ~ min(.x - 3 * .pggf.bw), 'max' = ~ max(.x + 3 * .pggf.bw), 'n' = length),
        .names = '.pggf.{.fn}'
      )
    ) |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of(c('facet_id', sample_cols)), !! pggf$facet_col, !! pggf$facet_row, !! pggf$facet_wrap)) |>
    dplyr::mutate(
      dplyr::across(
        {{extra_cols}},
        ~ if (length(unique(.x)) == 1) {
          unique(.x)
        } else {
          cli::cli_warn('Non-unique values in {.arg extra_cols}')
          paste0(.x, collapse = ',')
        }
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::nest_by(dplyr::pick(tidyselect::any_of(c('facet_id', sample_cols)), !! pggf$facet_col, !! pggf$facet_row, !! pggf$facet_wrap, {{extra_cols}})) |>
    dplyr::reframe(
      density = ifelse(
        !! tails == 'extend',
        list(density(data[[value_col]], from = data$.pggf.min[1], to = data$.pggf.max[1], n = {{samples}})),
        list(density(data[[value_col]], n = {{samples}}))
      ),
      !! value_col := density$x,
      density = density$y,
      n = length(data[[value_col]]),
      total_n = data$.pggf.n[1]
    ) |>
    dplyr::ungroup()
  
  # scale densities
  if (subdensity) {
    # normally, the density of each sample has an area of 1
    # with `subdensity = TRUE`, the sum of the area of all sample densities is 1
    # i.e. the subdensity takes into account the relative proportion of a sample in the total data
    ridgeline_data <- ridgeline_data |>
      dplyr::mutate(
        density = density * n / total_n
      )
  }
  
  # scale densities for non-numeric data
  if (! pggf[[paste0('numeric_', sample_cols[1])]]) {
    if (pggf$scales %in% c('free', paste0('free_', sample_cols[1]))) {
      ridgeline_data <- ridgeline_data |>
        dplyr::group_by(!! pggf[[paste0('facet_', ifelse(value_col == 'x', 'row', 'col'))]])
    }
    
    ridgeline_data <- ridgeline_data |>
      dplyr::mutate(
        density = density / max(density)
      ) |>
      dplyr::ungroup()
  }
  
  # select relevant columns and order data
  ridgeline_data <- ridgeline_data |>
    dplyr::select(tidyselect::any_of(c('facet_id', sample_cols)), !! pggf$facet_col, !! pggf$facet_row, !! pggf$facet_wrap, {{extra_cols}}, x, y, density) |>
    dplyr::arrange(
      dplyr::pick(tidyselect::any_of('facet_id')),
      dplyr::desc(dplyr::pick(tidyselect::any_of(sample_cols)))
    )
  
  # export data
  ridgeline_file <- ridgeline_data |>
    dplyr::select(tidyselect::any_of(c('facet_id', sample_cols)), {{extra_cols}}, x, y, !! paste0('density_', sample_cols[1]) := density) |>
    dplyr::mutate(
      file = paste0(pggf$filename, '_ridgeline.tsv')
    ) |>
    dplyr::nest_by(file) |>
    dplyr::summarise(
      content = readr::format_tsv(data),
      .groups = 'drop'
    ) |>
    tibble::deframe()
  
  # add file to pggf object
  pggf$files <- c(pggf$files, ridgeline_file)
  
  # replace original data with density data
  pggf$data <- ridgeline_data
  
  # add to ggplot
  if (! is.null(pggf$plot)) {
    if (orientation == 'x') {
      pggf$plot <- append(
        pggf$plot,
        ggplot2::geom_ribbon(
          mapping = ggplot2::aes(
            ymin = as.numeric(.data[[ifelse(pggf$numeric_y, 'y', 'y.orig')]]),
            ymax = as.numeric(.data[[ifelse(pggf$numeric_y, 'y', 'y.orig')]]) + density,
            fill = y,
            group = ordered(y, levels = sort(unique(y), decreasing = TRUE))
          ),
          orientation = value_col,
          color = 'black'
        )
      )
    } else {
      pggf$plot <- append(
        pggf$plot,
        ggplot2::geom_ribbon(
          mapping = ggplot2::aes(
            xmin = as.numeric(.data[[ifelse(pggf$numeric_x, 'x', 'x.orig')]]),
            xmax = as.numeric(.data[[ifelse(pggf$numeric_x, 'x', 'x.orig')]]) + density,
            fill = x,
            group = ordered(x, levels = sort(unique(x), decreasing = TRUE))
          ),
          orientation = value_col,
          color = 'black'
        )
      )
    }
  }
  
  # add to code
  if (! is.null(pggf$code)) {
    pggf$code <- append(
      pggf$code,
      '\t\\pggf_ridgeline()\n',
      after = length(pggf$code) - 1
    )
  }
  
  # call the new `pggf` object
  pggf
}


# logo plots
pggf_logo <- function(pggf, letters, align = NULL) {
  # make sure `pggf` is a `pggf` class object
  .assert_pggf(pggf, dims = 2)
  
  # make sure `x` and `y` are numeric
  if (! pggf$numeric_x || ! pggf$numeric_y) {
    cli::cli_abort('{.fn pggf_logo} only works with numeric x (position) and y (score) data')
  }
  # make sure `x` consists only of mathematical integers (does not actually have to be of integer type)
  if (any(pggf$data$x %% 1 != 0, na.rm = TRUE)) {
    cli::cli_abort('Only mathematical integers are allowed for {.arg x} with {.fn pggf_logo}')
  }
  
  # check arguments
  if (missing(letters)) {
    cli::cli_abort('Required argument {.arg letters} is missing')
  }
  if (! missing(align) && ! .pggf_is_null(substitute(align)) && ! align %in% c('left', 'center', 'right')) {
    cli::cli_abort('{.arg align} must be one of: "left", "center", or "right"')
  }
  
  # adjust x to align logos
  if (! missing(align) && ! .pggf_is_null(substitute(align))) {
    logo_data <- pggf$data |>
      dplyr::mutate(
        .pggf.line = seq_len(dplyr::n())
      ) |>
      tidyr::drop_na(x) |>
      dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
      dplyr::mutate(
        x.orig = x,
        x = x - min(x) + 1
      ) |>
      dplyr::ungroup()
    if (align == 'right') {
      logo_data <- logo_data |>
        dplyr::mutate(
          .pggf.max = max(x)
        ) |>
        dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
        dplyr::mutate(
          x = x + .pggf.max - max(x)
        ) |>
        dplyr::ungroup() |>
        dplyr::select(-.pggf.max)
    } else if (align == 'center') {
      logo_data <- logo_data |>
        dplyr::mutate(
          .pggf.max = max(x)
        ) |>
        dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
        dplyr::mutate(
          x = x + (.pggf.max - max(x)) %/% 2
        ) |>
        dplyr::ungroup() |>
        dplyr::select(-.pggf.max)
    }
    
    # update `pggf$data` with new x values
    pggf$data$x <- tibble::deframe(dplyr::select(logo_data, .pggf.line, x))[as.character(seq_along(pggf$data$x))]
    
    logo_data <- logo_data |>
      dplyr::select(-.pggf.line)
  } else {
    logo_data <- pggf$data |>
      tidyr::drop_na(x)
  }
  
  # prepare data
  logo_data <- logo_data |>
    dplyr::rename('letter' = {{letters}})
  
  # detect alphabet
  if (all(logo_data$letter %in% c(NA, 'A', 'C', 'G', 'T'))) {
    alphabet <- 'DNA'
  } else if (all(logo_data$letter %in% c(NA, 'A', 'C', 'G', 'U'))) {
    alphabet <- 'RNA'
  } else if (all(logo_data$letter %in% c(NA, strsplit('ACDEFGHIKLMNPQRSTVWY', '')[[1]]))) {
    alphabet <- 'protein'
  } else {
    alphabet <- 'rainbow'
  }
  
  # reshape data
  logo_data <- logo_data |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
    tidyr::complete(x, letter, fill = list(y = 0)) |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id')), x) |>
    dplyr::mutate(
      dplyr::across(tidyselect::any_of('x.orig'), ~ na.omit(.x)[1])
    ) |>
    dplyr::arrange(abs(y), .by_group = TRUE) |>
    dplyr::mutate(
      id = seq_len(dplyr::n())
    ) |>
    dplyr::ungroup() |>
    dplyr::select(tidyselect::any_of(c('facet_id', 'x', 'x.orig')), y, letter, id) |>
    tidyr::pivot_wider(
      names_from = id,
      values_from = c(y, letter)
    )
  
  # export data
  logo_file <- logo_data |>
    dplyr::mutate(
      file = paste0(pggf$filename, '_logo.tsv')
    ) |>
    dplyr::nest_by(file) |>
    dplyr::summarise(
      content = readr::format_tsv(data),
      .groups = 'drop'
    ) |>
    tibble::deframe()
  
  # add file to pggf object
  pggf$files <- c(pggf$files, logo_file)
  
  # save stack height and x limits for axis limit calculations
  logo_stack_range <- pggf$data |>
    dplyr::group_by(
      dplyr::pick(
        tidyselect::any_of(c('facet_id')),
        !! pggf$facet_col,
        !! pggf$facet_row,
        !! pggf$facet_wrap,
        x
      )
    ) |>
    dplyr::summarise(
      .pggf.min.y = sum(y[y < 0]),
      .pggf.max.y = sum(y[y > 0]),
      .groups = 'drop'
    ) |>
    dplyr::group_by(
      dplyr::pick(
        tidyselect::any_of(c('facet_id')),
        !! pggf$facet_col,
        !! pggf$facet_row,
        !! pggf$facet_wrap
      )
    ) |>
    dplyr::summarise(
      .pggf.min.x = min(x) - 0.5,
      .pggf.max.x = max(x) + 0.5,
      .pggf.min.y = max(.pggf.min.y),
      .pggf.max.y = max(.pggf.max.y),
      .groups = 'drop'
    ) |>
    tidyr::pivot_longer(
      tidyselect::starts_with(c('.pggf.min.', '.pggf.max.')),
      names_to = c('.pggf.name', '.value'),
      names_pattern = '\\.pggf\\.(m..)\\.(.)'
    ) |>
    dplyr::select(-.pggf.name)
  
  pggf$extra_coords <- pggf$extra_coords |>
    dplyr::bind_rows(
      logo_stack_range
    )
  
  # save alphabet length in axes file
  pggf$axes_extra <- pggf$data |>
    dplyr::group_by(dplyr::pick(tidyselect::any_of('facet_id'))) |>
    dplyr::summarise(
      `alphabet length` = dplyr::n_distinct({{letters}})
    )
  
  # add to ggplot
  if (! is.null(pggf$plot)) {
    cli::cli_warn(
      c(
        'Preview plot not implemented (yet) for {.fn pggf_logo}',
        'i' = 'Use {.arg preview_plot = FALSE} to silence this warning'
      )
    )
    pggf$plot <- NULL
  }
  
  # add to code
  if (! is.null(pggf$code)) {
    pggf$code <- append(
      pggf$code,
      paste0('\t\\pggf_logo(letter colors = ', alphabet, ')\n'),
      after = length(pggf$code) - 1
    )
  }
  
  # call the new `pggf` object
  pggf
}
