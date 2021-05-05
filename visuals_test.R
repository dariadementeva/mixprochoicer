# Visuals
plot_designD <- function(design){
  plot = ggtern(design,aes(V1,V2,V3))+
    geom_point(aes(V1,V2,V3), colour = "deepskyblue", size=10)+
    theme_showarrows()+
    labs(title = "D-optimal MPV Experiment",
      subtitle = "3 Ingredients, 30 runs, 2 PV",
      x = "Q1", xarrow = "Ingredient 1",
      y = "Q2", yarrow = "Ingredient 2",
      z = "Q3", zarrow = "Ingredient 3") +
    theme_latex()
  return(plot)
}

plot_designI <- function(design){
  plot = ggtern(design,aes(V1,V2,V3))+
    geom_point(aes(V1,V2,V3), colour = "deepskyblue", size=10)+
    theme_showarrows()+
    labs( title = "I-optimal Mixture Experiment",
      subtitle = "3 Ingredients, 30 runs, Order 3",
      x = "Q1", xarrow = "Ingredient 1",
      y = "Q2", yarrow = "Ingredient 2",
      z = "Q3", zarrow = "Ingredient 3") +
    theme_latex()
  return(plot)
}

# Gaussian Mixure Designs

# Order 1

random_design = create_random_initial_design(3,30,100000000)

compute_cox_direction(random_design[2,],3,30)

as.data.frame(compute_cox_direction(random_design[2,],3,100)) %>%
  ggtern() +
  geom_point(aes(V1,V2,V3), colour = "deepskyblue", size=0.2)+
  theme_showarrows()+
  labs(title = "Cox-Effect Direction: 100 Points",
    subtitle = "3 Ingredients, 30 runs",
    x = "Q1", xarrow = "Ingredient 1",
    y = "Q2", yarrow = "Ingredient 2",
    z = "Q3", zarrow = "Ingredient 3") +
  theme_latex()

optimized_designD = as.data.frame(coordinate_exchange_general(random_design,
  order=1,
  n_points=500,
  optimality_criterion="D",
  n_iter=100))

p1 <- plot_designD(optimized_designD)

optimized_designI = as.data.frame(coordinate_exchange_general(random_design,
  order=1,
  n_points=500,
  optimality_criterion="I",
  n_iter=100))

p2 <- plot_designI(optimized_designI)



# Order 2
optimized_designD_2 = as.data.frame(coordinate_exchange_general(random_design,
  order=2,
  n_points=500,
  optimality_criterion="D",
  n_iter=100))

p1_2 <- plot_designD(optimized_designD_2)

optimized_designI_2 = as.data.frame(coordinate_exchange_general(random_design,
  order=2,
  n_points=500,
  optimality_criterion="I",
  n_iter=100))

p2_2 <- plot_designI(optimized_designI_2)


# Order 3
#

optimized_designD_3 = as.data.frame(coordinate_exchange_general(random_design,
  order=3,
  n_points=100,
  optimality_criterion="D",
  n_iter=100))

p1_3 <- plot_designD(optimized_designD_3)

optimized_designI_3 = as.data.frame(coordinate_exchange_general(random_design,
  order=3,
  n_points=100,
  optimality_criterion="I",
  n_iter=100))

p2_3 <- plot_designI(optimized_designI_3)

grid.arrange(p1, p2, p1_2, p2_2, p1_3, p2_3, nrow = 3, ncol=2)


# Mixture-Choice Designs
design_MNL = create_random_initial_design_discrete(3,3,10,100000000)
q_MNL = 3
beta_MNL = rep(0, q - 1)
beta_MNL2 = rep(0, q - 1 + q * (q - 1)/2)
beta_MNL3 = rep(0,6)

Dopt_MNL_design= MNL_coordinate_exchange_d_optimal(design = design_MNL,
  order = 1,
  beta = beta_MNL,
  n_points = 100,
  n_iter = 100)

Dopt_MNL_design2= MNL_coordinate_exchange_d_optimal(design = design_MNL,
  order = 2,
  beta = beta_MNL2,
  n_points = 100,
  n_iter = 100)

Dopt_MNL_design3 = MNL_coordinate_exchange_d_optimal(design = design_MNL,
  order = 3,
  beta = beta_MNL3,
  n_points = 100,
  n_iter = 100)


mnl_design_array_to_dataframe = function(des_array, names = NULL){
  dim_X = dim(des_array)
  k = dim_X[1]
  S = dim_X[3]

  if(is.null(names)) names = c(paste0("c", 1:k), "choice_set")

  X_final_tbl = lapply(1:S, function(s){
    t(des_array[,,s]) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(choice_set = as.character(s))
  }) %>%
    dplyr::bind_rows() %>%
    purrr::set_names(names)

  return(X_final_tbl)

}

mnl1 = mnl_design_array_to_dataframe(Dopt_MNL_design, names = c("V1", "V2", "V3", "choice_set"))
mnl1_g = ggtern(mnl1, aes(V1,V2,V3), colour = mnl1$choice_set)+
  geom_point(aes(V1,V2,V3), colour = "deepskyblue", size=3)+
  theme_showarrows()+
  labs(x = "Q1", xarrow = "Ingredient 1",
    y = "Q2", yarrow = "Ingredient 2",
    z = "Q3", zarrow = "Ingredient 3") +
  theme_nomask()
  theme_latex()

mnl2 = mnl_design_array_to_dataframe(Dopt_MNL_design2, names = c("V1", "V2", "V3", "choice_set"))
mnl2_g = ggtern(mnl2, aes(V1,V2,V3), colour = mnl1$choice_set)+
  geom_point(aes(V1,V2,V3), colour = "deepskyblue", size=3)+
  theme_showarrows()+
  labs(x = "Q1", xarrow = "Ingredient 1",
    y = "Q2", yarrow = "Ingredient 2",
    z = "Q3", zarrow = "Ingredient 3") +
  theme_nomask()
  theme_latex()


mnl3 = mnl_design_array_to_dataframe(ccc)
mnl3_g =ggtern(mnl3, aes(c1,c2, c3), colour = mnl1$choice_set)+
  geom_point(aes(c1,c2,c3), colour = "deepskyblue", size=3)+
  theme_showarrows()+
  labs(x = "Q1", xarrow = "Ingredient 1",
    y = "Q2", yarrow = "Ingredient 2",
    z = "Q3", zarrow = "Ingredient 3") +
  theme_latex()+
  theme_nomask()

grid.arrange(mnl1_g, mnl2_g, mnl3_g, nrow = 3, ncol=1, top = "D-optimal Mixture Choice Experiment:
  3 Ingredients, 3 Alternatives, 10 Choice Sets, Order 1, 2, 3")


MPV_design = MPV_random_design(3,30,2)
optimized_design = MPV_coordinate_exchange_d_optimal(MPV_design,3,2,100,100)


# ------------- COMPARISON ------------------------ #

# 1

mm1 = MNL_MPV_create_moment_matrix(3,1,4)
all.equal(mm1, mnl_moments_matrix_q3_npv1_order4_pvboundsminus1to1) #TRUE

# 2

mm2 = MNL_MPV_create_moment_matrix(3,2,4)
all.equal(mm2, mnl_moments_matrix_q3_npv2_order4_pvboundsminus1to1) #TRUE

# 3

mm3 = MNL_MPV_create_moment_matrix(3,3,4)
all.equal(mm3, mnl_moments_matrix_q3_npv3_order4_pvboundsminus1to1) #TRUE

# 4

mm4 = MNL_MPV_create_moment_matrix(3,4,4)
all.equal(mm4, mnl_moments_matrix_q3_npv4_order4_pvboundsminus1to1) #TRUE

# 5

mm5 = MNL_MPV_create_moment_matrix(3,5,4)
all.equal(mm5,mnl_moments_matrix_q3_npv5_order4_pvboundsminus1to1) #TRUE

# 6

mm5 = MNL_MPV_create_moment_matrix(3,5,4)
all.equal(mm5,mnl_moments_matrix_q3_npv5_order4_pvboundsminus1to1) #TRUE

# 6

mm6 = MNL_MPV_create_moment_matrix(3,0,1)
all.equal(mm6,mnl_moments_matrix_q3_order1_npv0) #TRUE

# 7

mm7 = MNL_MPV_create_moment_matrix(3,0,2)
all.equal(mm7,mnl_moments_matrix_q3_order2_npv0) #TRUE

# 8

mm8 = MNL_MPV_create_moment_matrix(3,0,3)
all.equal(mm8,mnl_moments_matrix_q3_order3_npv0) #TRUE

# 9

mm9 = MNL_MPV_create_moment_matrix(4,1,4)
all.equal(mm9,mnl_moments_matrix_q4_npv1_order4_pvboundsminus1to1) #TRUE

# 10

mm10 = MNL_MPV_create_moment_matrix(4,2,4)
all.equal(mm10,mnl_moments_matrix_q4_npv2_order4_pvboundsminus1to1) #TRUE

# 11

mm11 = MNL_MPV_create_moment_matrix(4,3,4)
all.equal(mm11,mnl_moments_matrix_q4_npv3_order4_pvboundsminus1to1) #TRUE


# 12

mm12 = MNL_MPV_create_moment_matrix(4,4,4)
all.equal(mm12,mnl_moments_matrix_q4_npv4_order4_pvboundsminus1to1) #TRUE

# 13

mm13 = MNL_MPV_create_moment_matrix(4,5,4)
all.equal(mm13,mnl_moments_matrix_q4_npv5_order4_pvboundsminus1to1) #TRUE

# 14

mm14 = MNL_MPV_create_moment_matrix(4,0,1)
all.equal(mm14,mnl_moments_matrix_q4_order1_npv0) #TRUE

# 15

mm15 = MNL_MPV_create_moment_matrix(4,0,2)
all.equal(mm15,mnl_moments_matrix_q4_order2_npv0) #TRUE

# 16

mm16 = MNL_MPV_create_moment_matrix(4,0,3)
all.equal(mm16,mnl_moments_matrix_q4_order3_npv0) #TRUE


# 17

mm17 = MNL_MPV_create_moment_matrix(5,1,4)
all.equal(mm17,mnl_moments_matrix_q5_npv1_order4_pvboundsminus1to1) #TRUE

# 18

mm18 = MNL_MPV_create_moment_matrix(5,2,4)
all.equal(mm18,mnl_moments_matrix_q5_npv2_order4_pvboundsminus1to1) #TRUE

# 19

mm19 = MNL_MPV_create_moment_matrix(5,3,4)
all.equal(mm19,mnl_moments_matrix_q5_npv3_order4_pvboundsminus1to1) #TRUE


# 20

mm20 = MNL_MPV_create_moment_matrix(5,4,4)
all.equal(mm20,mnl_moments_matrix_q5_npv4_order4_pvboundsminus1to1) #TRUE

# 21

mm21 = MNL_MPV_create_moment_matrix(5,5,4)
all.equal(mm21,mnl_moments_matrix_q5_npv5_order4_pvboundsminus1to1) #TRUE


# 22

mm22 = MNL_MPV_create_moment_matrix(5,0,1)
all.equal(mm22,mnl_moments_matrix_q5_order1_npv0) #TRUE

# 23

mm23 = MNL_MPV_create_moment_matrix(5,0,2)
all.equal(mm23,mnl_moments_matrix_q5_order2_npv0) #TRUE

# 24

mm24 = MNL_MPV_create_moment_matrix(5,0,3)
all.equal(mm24,mnl_moments_matrix_q5_order3_npv0) #TRUE















