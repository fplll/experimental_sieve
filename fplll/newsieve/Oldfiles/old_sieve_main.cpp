// Code from old sieve_main.cpp

#if 0
  dot_time = 0;
  dot_num = 0;
  count_bad = 0;
#endif

  alg = 2;

  /* parse */
  if (argc == 1)
  {
    main_usage(argv[0]);
    return -1;
  }
  while ((option = getopt(argc, argv, "a:f:r:t:s:b:v")) != -1)
  {
    switch (option)
    {
    case 'a':
      alg = atoi(optarg);
      if (alg != 2 && alg != 3 && alg != 4)
      {
        cout << " Error, only support 2-, 3- and 4-sieve" << endl;
        exit(1);
      }
    case 'f':
      input_file_name = optarg;
      flag_file       = true;
      break;
    case 'r':
      dim       = atoi(optarg);
      flag_file = false;
      break;
    case 's':
      seed = atoi(optarg);
      break;
    case 'b':
      bs = atoi(optarg);
      break;
    case 'v':
      flag_verbose = true;
      break;
    case 't':
      // ngoal_norm = atol(optarg);
      cout << optarg << endl;
      goal_norm_s = optarg;
      break;
    case 'h':
      main_usage(argv[0]);
      return -1;
    case '?':
      main_usage(argv[0]);
      return -1;
    case ':':
      main_usage(argv[0]);
      return -1;
    }
  }

  /* set lattice */
  ZZ_mat<mpz_t> B;
  if (flag_file)
  {
    ifstream input_file(input_file_name);
    if (input_file.is_open())
    {
      input_file >> B;
      input_file.close();
    }
    else
    {
      cin >> B;
    }
    if (flag_verbose)
    {
      cout << "# [info] reading lattice of dimension " << B.get_rows() << "x" << B.get_cols()
           << endl;
    }
  }
  else
  {
    if (flag_verbose)
    {
      cout << "# [info] generating random lattice of dimension " << dim << endl;
    }
    srand(time(NULL));
    B.resize(dim, dim);
    //B.gen_trg(1.1);
    BTest.gen_qary_prime(1, 10*dim);
  }

  /* set targeted norm */
  Z_NR<mpz_t> goal_norm, max;
  if (goal_norm_s != NULL)
  {
    goal_norm.set_str(goal_norm_s);
  }
  if (goal_norm < 0)
    goal_norm = 0;
  if (flag_verbose)
    cout << "# [info] goal norm^2 is " << goal_norm << endl;

  /* preprocessing of basis */
  clock_t stime = clock();
  if (bs > 0)
    bkz_reduction(B, bs, BKZ_DEFAULT, FT_DEFAULT, 0);
  else
      lll_reduction(B, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER);

  clock_t etime = clock();
  double secs   = (etime - stime) / (double)CLOCKS_PER_SEC;
  if (flag_verbose)
  {
    if (bs > 0)
      cout << "# [info] BKZ took time " << secs << " s" << endl;
    else
      cout << "# [info] LLL took time " << secs << " s" << endl;
  }
  // cout << B << endl;

  /* decide integer type */
  stime = clock();
  max   = B.get_max();

#if 1
  if (max < std::numeric_limits<int>::max())
  {
    long goal_norm_l = abs(goal_norm.get_si());
    Z_NR<long> goal_norm_lt;
    goal_norm_lt = goal_norm_l;
    ZZ_mat<long> B2(B.get_rows(), B.get_cols());
    for (int i = 0; i < B.get_rows(); i++)
      for (int j = 0; j < B.get_cols(); j++)
        B2(i, j) = B(i, j).get_si();
    main_run_sieve<long>(B2, goal_norm_lt, alg, flag_verbose, seed);
  }
  else
#endif
 //   main_run_sieve<mpz_t>(B, goal_norm, alg, flag_verbose, seed); //commented out for test

  etime = clock();
  secs  = (etime - stime) / (double)CLOCKS_PER_SEC;
  if (flag_verbose)
  {
    cout << "# [info] sieve took time " << secs << " s" << endl;
/* dot product time */
#if 0
    cout << "# [info] dot_time " << dot_time << endl;
    cout << "# [info] dot_num " << dot_num << endl;
    cout << "# [info] dot_time/dot_number " << (double) dot_time/dot_num << endl;
  }
#endif
