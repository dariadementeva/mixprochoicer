# Fish Patty Experiment

# D-optimality NL_MPV_random_design(q = 3, J = 2, S = 20, m = 1, m_bounds =
# c(-1,1), seed = 100)

set.seed(5081997)

fish_patty_random_starting_design = MNL_MPV_random_design(q = 3, J = 2, S = 28, m = 3, 
    m_bounds = c(-1, 1), seed = 5081997)

starting_fish_patty_dopt_value = MNL_MPV_compute_d_optimality(design = fish_patty_random_starting_design, 
    beta = rep(0, 20), m = 3, order = 4)  # 60.02469

d_optimized_fish_patty_design = MNL_MPV_coordinate_exchange_general(design = fish_patty_random_starting_design, 
    m = 3, order = 4, beta = rep(0, 20), n_points = 500, n_iter = 50, optimality_criterion = "D", 
    m_bounds = c(1, 1))

optimized_fish_patty_dopt_value = MNL_MPV_compute_d_optimality(design = d_optimized_fish_patty_design, 
    beta = rep(0, 20), m = 3, order = 4)  # -29.65346

# Checking whether the optimality value has been minimized
starting_fish_patty_dopt_value > optimized_fish_patty_dopt_value  #TRUE

# Fish Patty Experiment

# I-optimality m = 3, order = 4, beta = rep(0, 20), m_bounds = c(1, 1))

starting_fish_patty_iopt_value = MNL_MPV_compute_i_optimality(design = fish_patty_random_starting_design, 
    m = 3, order = 4, beta = rep(0, 20), m_bounds = c(1, 1))  #223.8802

i_optimized_fish_patty_design = MNL_MPV_coordinate_exchange_general(design = fish_patty_random_starting_design, 
    m = 3, order = 4, beta = rep(0, 20), n_points = 1000, n_iter = 10, optimality_criterion = "I", 
    m_bounds = c(1, 1))

optimized_fish_patty_iopt_value = MNL_MPV_compute_d_optimality(design = i_optimized_fish_patty_design, 
    m = 3, order = 4, beta = rep(0, 20), m_bounds = c(1, 1))

# Checking whether the optimality value has been minimized
starting_fish_patty_iopt_value > optimized_fish_patty_iopt_value


# PARALLEL PROCESSING, D-OPT

d_opt_list = mclapply(1:18, function(i) {

    print(i)

    fish_patty_random_starting_design_1 = MNL_MPV_random_design(q = 3, J = 2, S = 28, 
        m = 3, m_bounds = c(-1, 1), seed = i)

    d_optimized_fish_patty_design_1 = MNL_MPV_coordinate_exchange_general(design = fish_patty_random_starting_design_1, 
        m = 3, order = 4, beta = rep(0, 20), n_points = 100, n_iter = 10, optimality_criterion = "D", 
        m_bounds = c(1, 1))

    optimized_fish_patty_dopt_value = MNL_MPV_compute_d_optimality(design = d_optimized_fish_patty_design_1, 
        beta = rep(0, 20), m = 3, order = 4)

    out_list = list(fish_patty_random_starting_design_1, d_optimized_fish_patty_design_1, 
        optimized_fish_patty_dopt_value)
    return(out_list)
}, mc.cores = 6)

d_optimal_design_for_fish_patty_experiment = d_opt_list[[15]]
d_optimal_design_for_fish_patty_experiment_initial_design_array = d_opt_list[[15]][1]
d_optimal_design_for_fish_patty_experiment_design_array = d_opt_list[[15]][2]
d_optimal_design_for_fish_patty_experiment_d_opt_value = d_opt_list[[15]][3]


# PARALLEL PROCESSING, I-OPT

i_opt_list = mclapply(1:18, function(i) {

    print(i)

    fish_patty_random_starting_design_2 = MNL_MPV_random_design(q = 3, J = 2, S = 28, 
        m = 3, m_bounds = c(-1, 1), seed = i)

    i_optimized_fish_patty_design_1 = MNL_MPV_coordinate_exchange_general(design = fish_patty_random_starting_design_2, 
        m = 3, order = 4, beta = rep(0, 20), n_points = 100, n_iter = 10, optimality_criterion = "I", 
        m_bounds = c(1, 1))

    optimized_fish_patty_iopt_value = MNL_MPV_compute_i_optimality(design = i_optimized_fish_patty_design_1, 
        m = 3, order = 4, beta = rep(0, 20), m_bounds = c(1, 1))

    out_list = list(fish_patty_random_starting_design_2, i_optimized_fish_patty_design_1, 
        optimized_fish_patty_iopt_value)
    return(out_list)
}, mc.cores = 6)


i_optimal_design_for_fish_patty_experiment = i_opt_list[[10]]
i_optimal_design_for_fish_patty_experiment_initial_design_array = i_opt_list[[10]][1]
i_optimal_design_for_fish_patty_experiment_design_array = i_opt_list[[10]][2]
i_optimal_design_for_fish_patty_experiment_i_opt_value = i_opt_list[[10]][3]


# VISUALS: D-Optimal

# convert from list to an array
FP_D_OPT_DES = array(as.numeric(unlist(d_optimal_design_for_fish_patty_experiment_design_array)), 
    dim = c(6, 2, 28))

# convert from an array to data frame
FP_D_OPT_DES_DF = mnl_design_array_to_dataframe(FP_D_OPT_DES)


FP_D_OPT_DES_DF <- FP_D_OPT_DES_DF %>%
    rename(`Deep Frying Time` = c4, `Cooking Temperature` = c5, `Cooking Time` = c6)


plotD = FP_D_OPT_DES_DF %>%
    pivot_longer("Deep Frying Time":"Cooking Time") %>%
    ggtern() + geom_point(aes(x = c1, y = c2, z = c3, color = value)) + scale_color_gradient(low = "hotpink", 
    high = "hotpink4") + facet_wrap(~name) + theme_light() + theme_nomask() + theme_showarrows() + 
    labs(title = "Fish Patty MPV Choice Experiment. D-optimality Value = -30.3244", 
        subtitle = "3 Ingredients, 3 Process Variables", x = "Q1", xarrow = "Mullet", 
        y = "Q2", yarrow = "Sheephead", z = "Q3", zarrow = "Croacker")

ggsave("plotD.png")

plotD_CS = FP_D_OPT_DES_DF %>%
    pivot_longer("Deep Frying Time":"Cooking Time") %>%
    ggtern() + geom_point(aes(x = c1, y = c2, z = c3, color = choice_set)) + # scale_color_gradient(low = 'hotpink', high = 'hotpink4')+
facet_wrap(~name) + theme_light() + theme_nomask() + theme_showarrows() + theme(legend.position = "bottom") + 
    labs(title = "Fish Patty MPV Choice Experiment. D-optimality Value = -30.3244", 
        subtitle = "3 Ingredients, 3 Process Variables", x = "Q1", xarrow = "Mullet", 
        y = "Q2", yarrow = "Sheephead", z = "Q3", zarrow = "Croacker")

ggsave("plotD_CS.png")



# VISUALS: D-Optimal EXTRA
d_optimal_design_for_fish_patty_experiment_design_array = d_opt_list[[18]][2]

# convert from list to an array
FP_D_OPT_DES = array(as.numeric(unlist(d_optimal_design_for_fish_patty_experiment_design_array)), 
    dim = c(6, 2, 28))

# convert from an array to data frame
FP_D_OPT_DES_DF = mnl_design_array_to_dataframe(FP_D_OPT_DES)


FP_D_OPT_DES_DF <- FP_D_OPT_DES_DF %>%
    rename(`Deep Frying Time` = c4, `Cooking Temperature` = c5, `Cooking Time` = c6)


plotD = FP_D_OPT_DES_DF %>%
    pivot_longer("Deep Frying Time":"Cooking Time") %>%
    ggtern() + geom_point(aes(x = c1, y = c2, z = c3, color = value)) + scale_color_gradient(low = "hotpink", 
    high = "hotpink4") + facet_wrap(~name) + theme_light() + theme_nomask() + theme_showarrows() + 
    labs(title = "Fish Patty MPV Choice Experiment. D-optimality Value = -30.3244", 
        subtitle = "3 Ingredients, 3 Process Variables", x = "Q1", xarrow = "Mullet", 
        y = "Q2", yarrow = "Sheephead", z = "Q3", zarrow = "Croacker")

ggsave("plotD.png")

plotD_CS = FP_D_OPT_DES_DF %>%
    pivot_longer("Deep Frying Time":"Cooking Time") %>%
    ggtern() + geom_point(aes(x = c1, y = c2, z = c3, color = choice_set)) + # scale_color_gradient(low = 'hotpink', high = 'hotpink4')+
facet_wrap(~name) + theme_light() + theme_nomask() + theme_showarrows() + theme(legend.position = "bottom") + 
    labs(title = "Fish Patty MPV Choice Experiment. D-optimality Value = -30.3244", 
        subtitle = "3 Ingredients, 3 Process Variables", x = "Q1", xarrow = "Mullet", 
        y = "Q2", yarrow = "Sheephead", z = "Q3", zarrow = "Croacker")

ggsave("plotD_CS.png")










# VISUALS: I-Optimal



# convert from list to an array
FP_I_OPT_DES = array(as.numeric(unlist(i_optimal_design_for_fish_patty_experiment_design_array)), 
    dim = c(6, 2, 28))

# convert from an array to data frame
FP_I_OPT_DES_DF = mnl_design_array_to_dataframe(FP_I_OPT_DES)

View(FP_I_OPT_DES_DF %>%
    group_by(choice_set) %>%
    summarise(count = n_distinct(FP_I_OPT_DES_DF$c1:FP_I_OPT_DES_DF$c3)))




FP_I_OPT_DES_DF <- FP_I_OPT_DES_DF %>%
    rename(`Deep Frying Time` = c4, `Cooking Temperature` = c5, `Cooking Time` = c6)


plotI = FP_I_OPT_DES_DF %>%
    pivot_longer("Deep Frying Time":"Cooking Time") %>%
    ggtern() + geom_point(aes(x = c1, y = c2, z = c3, color = value)) + scale_color_gradient(low = "hotpink", 
    high = "hotpink4") + facet_wrap(~name) + theme_light() + theme_nomask() + theme_showarrows() + 
    labs(title = "Fish Patty MPV Choice Experiment. I-optimality Value = 2.81", subtitle = "3 Ingredients, 3 Process Variables", 
        x = "Q1", xarrow = "Mullet", y = "Q2", yarrow = "Sheephead", z = "Q3", zarrow = "Croacker")

ggsave("plotI.png")

plotI_CS = FP_I_OPT_DES_DF %>%
    pivot_longer("Deep Frying Time":"Cooking Time") %>%
    ggtern() + geom_point(aes(x = c1, y = c2, z = c3, color = choice_set)) + # scale_color_gradient(low = 'hotpink', high = 'hotpink4')+
facet_wrap(~name) + theme_light() + theme_nomask() + theme_showarrows() + theme(legend.position = "bottom") + 
    labs(title = "Fish Patty MPV Choice Experiment. I-optimality Value = 2.7411", 
        subtitle = "3 Ingredients, 3 Process Variables", x = "Q1", xarrow = "Mullet", 
        y = "Q2", yarrow = "Sheephead", z = "Q3", zarrow = "Croacker")

ggsave("plotI_CS.png")





# EXTRA I


i_optimal_design_FP_extra = i_opt_list[[18]][2]

# convert from list to an array
FP_I_OPT_DES_extra = array(as.numeric(unlist(i_optimal_design_FP_extra)), dim = c(6, 
    2, 28))

# convert from an array to data frame
FP_I_OPT_DES_DF_extra = mnl_design_array_to_dataframe(FP_I_OPT_DES_extra)


FP_I_OPT_DES_DF_extra <- FP_I_OPT_DES_DF_extra %>%
    rename(`Deep Frying Time` = c4, `Cooking Temperature` = c5, `Cooking Time` = c6)


plotI_extra = FP_I_OPT_DES_DF_extra %>%
    pivot_longer("Deep Frying Time":"Cooking Time") %>%
    ggtern() + geom_point(aes(x = c1, y = c2, z = c3, color = value)) + scale_color_gradient(low = "hotpink", 
    high = "hotpink4") + facet_wrap(~name) + theme_light() + theme_nomask() + theme_showarrows() + 
    labs(title = "Fish Patty MPV Choice Experiment. I-optimality Value = 2.81", subtitle = "3 Ingredients, 3 Process Variables", 
        x = "Q1", xarrow = "Mullet", y = "Q2", yarrow = "Sheephead", z = "Q3", zarrow = "Croacker")

ggsave("plotI_extra.png")

















