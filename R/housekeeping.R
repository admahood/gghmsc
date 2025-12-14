# housekeeping script
ignore_unused_imports <- function(){
  utils::globalVariables

}

# solving the note on global variables (a symptom of using dplyr a lot)
utils::globalVariables(c('Support', 'env_var', 'Mean', 'value', 'prevalence',
                         'Species', 'species', 'var1', 'order_f', 'Species_f',
                         'var', 'sp', 'Chain', 'Parameter', 'prev_pct',
                         'Iteration', 'fg', 'origin', 'x', 'p_equiv',
                         'median_value', 'xintercept', 'Point est.', '.',
                         'Trait', 'variable', 'variable_pct', 'name',
                         'sp_f', 'trait', 'value_pct', 'x1', 'x2'))
