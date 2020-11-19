test_that("Walk-through (v16.4) of the model code", {
  # skip("temp skip")

  rm(list = ls())
  source(paste0(getwd(), "/common.R"), local = environment())
  init(e = environment())

  # file_path <- paste0(getwd(), "/data/templates_v16.2/Template_CoMoCOVID-19App_r_v_cpp.xlsx")
  # file_path <- paste0(getwd(), "/data/templates_v16.3/Testtemplate3.xlsx")
  file_path <- paste0(getwd(), "/data/templates_v16.4/Template_CoMoCOVID-19App_v17.xlsx")

  if (!exists("inputs", mode = "function")) {
    # source(paste0(getwd(), CORE_FILE), local = environment())
    # print("before source")
    # source(paste0(getwd(), "/v16.3.core.R"), local = environment())
    source(paste0(getwd(), "/v16.4.core.mod.R"), local = environment())
    # print("after source")
  }

  # stop()

  check_parameters_list_for_na(parameters)

  # environment(check_mortality_count) <- environment()

  # p_value_list <- seq(0.0, 0.1, by = 0.025)
  p_value_list <- seq(0.1, 0.1)

  scenario_list <- list(
    # vectors0, # Baseline
    vectors0 # Baseline
    # vectors   # Hypothetical
  )

  for (ss in scenario_list) {
    for (pp in p_value_list) {
      parameters["p"] <- pp

      param_vector <- parameters
      # param_vector[parameters_noise] <- parameters[parameters_noise]
      #   + rnorm(
      #       length(parameters_noise),
      #       mean = 0,
      #       sd = noise * abs(parameters[parameters_noise])
      #     )

      # start_time <- Sys.time()
      covidOdeCpp_reset()
      output_message <- capture_output(
        out_cpp <- ode(
                    y = Y, times = times, method = "euler", hini = 0.05,
                    func = covidOdeCpp, parms = param_vector,
                    input = ss, A = A,
                    contact_home = contact_home,
                    contact_school = contact_school,
                    contact_work = contact_work,
                    contact_other = contact_other,
                    popbirth_col2 = popbirth[, 2],
                    popstruc_col2 = popstruc[, 2],
                    ageing = ageing,
                    ifr_col2 = ifr[, 2],
                    ihr_col2 = ihr[, 2],
                    mort_col = mort,
                    age_group_vectors = age_group_vectors
                    )

      )
      # elapsed_time <- Sys.time() - start_time

      # print("Rcpp version time:")
      # print(elapsed_time)

      expect_equal(output_message, "covidOdeCpp: splinefuns updated")

      processed_cpp_results <- process_ode_outcome_mortality(out_cpp, ss, param_vector)

      expect_equal(
        processed_cpp_results$total_reportable_deaths_end,
        processed_cpp_results$total_cm_deaths_end,
        tolerance = 0.1,
        scale = processed_cpp_results$total_cm_deaths_end
      )


      # start_time <- Sys.time()
      expect_silent(
        out_r <- ode(
                  y = Y, times = times, method = "euler", hini = 0.05,
                  func = covid, parms = param_vector, input = ss
                  )
      )
      # elapsed_time <- Sys.time() - start_time

      # print("R version time:")
      # print(elapsed_time)

      processed_r_results <- process_ode_outcome_mortality(out_r, ss, param_vector)

      expect_equal(
        processed_r_results$total_reportable_deaths_end,
        processed_r_results$total_cm_deaths_end,
        tolerance = 0.1,
        scale = processed_r_results$total_cm_deaths_end
      )

      match_processed_outputs(
          output_a = processed_cpp_results,
          output_b = processed_r_results,
          tlr = 0.0001
      )

      # # sss = 1
      # # write.csv(out_cpp, paste0("out_cpp_",sss,"_",parameters["p"],".csv"),row.names = FALSE)
      # # write.csv(out_r, paste0("out_r_",sss,"_",parameters["p"],".csv"),row.names = FALSE)

      match_outputs(
        output_a = out_r,
        output_b = out_cpp,
        tlr = 0.0001,
        smp = 1000
      )

    }
  }

})