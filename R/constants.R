ra_drugs = c('hydroxychloroquine', 'sulfasalazine', 'methotrexate', 'leflunomide',
						 'adalimumab', 'infliximab', 'etanercept', 'certolizumab', 'golimumab',
						 'abatacept',
						 'tocilizumab', 'sarilumab',
						 'anakinra',
						 'rituximab',
						 'tofacitinib', 'baricitinib', 'upadacitinib', 'filgotinib')

ra_drug_regex = paste0('(', paste(ra_drugs, collapse = ')|('), ')')
