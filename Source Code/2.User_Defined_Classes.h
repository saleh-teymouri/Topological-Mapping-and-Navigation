class Landmark: public AStar::Node<Landmark,double> {
public:

	Point coord;
	int id;
	unordered_set<Landmark*> neighbors;
	unordered_map<Landmark*, int> neighbors_occlusion_count;
	double distances;
	int obs_count;
	int g_score;
	int rbt_id;
	bool NearOBS = false;

	Landmark () { }

	Landmark (int x, int y) {
		coord = Point(x, y);
	}

	Landmark (int x, int y, int i) {
		coord = Point(x, y);
		id = i;
	}

	bool operator==(const Landmark & n) const { return (id == n.id); }

	int getHashBin (void) { return (abs(((int)id>>4) + ((int)id<<3))); }
};

bool Landmark_Comparator (const Landmark & l1, Landmark & l2) {
return (l1.distances < l2.distances); }

class dpoint {
public:

	Point coord;
	double angle;
	double sensor_radius;
	double sensor_angle;
	int rbt_id;
	vector<Landmark> observing_landmarks;
	vector<int> nearest_landmark_ids;
	vector<int> gscore_to_clusters;
	bool Assigned = false;

	dpoint (double x, double y, double a) {
		coord = Point(x, y);
		angle = a;
	}

	dpoint (double x, double y, double a, double sr, double sa, int id) {
		coord = Point(x, y);
		angle = a;		
		sensor_radius = sr;
		sensor_angle = sa;
		rbt_id = id;
	}
};

class edges {
public:

	double homology_value;
	int edge_id;

	edges () { }

	edges (double h, int id) {
		homology_value = h;
		edge_id = id;
	}
};

bool edges_Comparator (const edges & e1, edges & e2) { 
return (e1.homology_value < e2.homology_value); }

/*************************************************************************************************/
/** Global Variables **//** Start **/

vector<Landmark> landmarks;
vector<dpoint> observations;
vector<edges> H;
unordered_set<int> cc_ids;
vector<vector<Landmark>> clusters;
vector<unordered_set<int>> clusters_ids;
vector<unordered_set<int>> obs_list(1);

/** Global Variables **//** End **/
/*************************************************************************************************/

int nCr (int n,int r); //3.User_Defined_Functions.h:28
int locate_Jump (vector<double> data); //3:535
void identify_clusters (); //3:616
void draw_vector (arma::Mat<double> x, Mat & img);//3:655
void assign_directions (arma::Mat<double> & x, unordered_set<int> edge_ids);//3:711

/*************************************************************************************************/

template <class v_type> class simplex {
public:

	unordered_set<v_type> S;
	int id;
	bool c2_count = false;
};

template <class T1> class Hasher {
public:

	size_t operator() (const simplex<T1>& obj) const {
		size_t sum = 0;
		for (auto n = obj.S.begin(); n != obj.S.end(); ++n) 
			sum += *n;
		return sum;
	}
};

template <class T2> class Comparator {
public:

	bool operator() (const simplex<T2>& obj1, const simplex<T2>& obj2) const {

		if (obj1.S.size() != obj2.S.size())
			return false;
		vector<T2> v1 (obj1.S.begin(), obj1.S.end());
		sort (v1.begin(), v1.end());
		vector<T2> v2 (obj2.S.begin(), obj2.S.end());
		sort (v2.begin(), v2.end());
		for(int n = 0; n < v1.size(); ++n) {
			if (v1.at(n) != v2.at(n))
				return false;
		}
		return true;
	}
};

template <class v_type> class simplicial_complex {
public:

	simplex<v_type> C0;
	vector<unordered_set<simplex<v_type>, Hasher<v_type>, Comparator<v_type>>> C;
	vector<simplex<v_type>> id_to_simplex;
	vector<simplex<v_type>> c2_pointers;

	void create_observation (simplex<v_type> observed_simplex) {

		if (observed_simplex.S.size() == 0)
			return;

		if (C.size() < observed_simplex.S.size())
			C.resize(observed_simplex.S.size());

		if (C[observed_simplex.S.size()-1].find(observed_simplex) ==
			C[observed_simplex.S.size()-1].end()) {

			if (observed_simplex.S.size() == 1) {
				for (auto it = observed_simplex.S.begin(); it != observed_simplex.S.end(); ++it) {
					if (*it > 0)
						C0.S.insert(*it);
				}

				return;
			}

			if (observed_simplex.S.size() == 2) {
				observed_simplex.id = C[1].size();
				id_to_simplex.push_back(observed_simplex);
			}

			if (observed_simplex.S.size() == 3) {
				observed_simplex.id = C[2].size();
				c2_pointers.push_back(observed_simplex);
/*				bool BRK1 = false;
				bool BRK2 = false;
				for (auto it = observed_simplex.S.begin(); it != observed_simplex.S.end(); ++it) {
					for (auto i = landmarks[*it].neighbors.begin();
					i != landmarks[*it].neighbors.end(); ++i) {
						simplex<v_type> observed_simplex_copy = observed_simplex;
						simplex<v_type> new_simplex;
						observed_simplex_copy.S.erase(*it);

						if (observed_simplex_copy.S.find((**i).id) !=
						observed_simplex_copy.S.end())
							continue;

						bool CNT = false;
						for (auto it2 = observed_simplex_copy.S.begin();
						it2 != observed_simplex_copy.S.end(); ++it2) {
							if (landmarks[*it2].neighbors.find(*i) ==
							landmarks[*it2].neighbors.end()) {
								CNT = true;
								break;
							}

							else {
								new_simplex.S.insert ((**i).id);
								new_simplex.S.insert (*it2);
							}
						}

						if (CNT)
							continue;

						if (C[2].find(new_simplex) != C[2].end()) {
							BRK1 = true;
							BRK2 = true;
							break;
						}

						if (BRK1)
							break;
					}

					if (BRK2)
						break;
				}

				if (!observed_simplex.c2_count && !BRK1) {
					observed_simplex.c2_count = true;
					C2_size++;
				}
*/			}

			C[observed_simplex.S.size()-1].insert(observed_simplex);
			sc_count++;

			simplex<v_type> temp;
			for (auto it = observed_simplex.S.begin(); it != observed_simplex.S.end(); ++it) {
				temp = observed_simplex;
				temp.S.erase(*it);
				create_observation (temp);
			}
		}

		else
			return;
	}

	void redundant_simplex (simplex<v_type> observed_simplex) {

		int sum_f(0);
		sum_f = nCr(observed_simplex.S.size(), 3) - c2_pointers.size();

		for (int it = 0; it < c2_pointers.size(); ++it) {
			if (c2_pointers[it].c2_count)
				sum_f++;
		}

		int c2_temp_size = (observed_simplex.S.size()-2) - sum_f;
		if (c2_temp_size > 0)
			C2_size += c2_temp_size;

		for (int it = 0; it < c2_pointers.size(); ++it) {
			C[2].erase(c2_pointers[it]);
			c2_pointers[it].c2_count = true;
			C[2].insert(c2_pointers[it]);
		}
	}

	void Homology () {

		arma::sp_mat D1(landmarks.size(), C[1].size());
		arma::sp_mat D2(C[1].size(), C[2].size());

//		arma::Mat<double> D1(landmarks.size(), C[1].size());
//		arma::Mat<double> D2(C[1].size(), C[2].size());
	
		for (auto it = C[1].begin(); it != C[1].end(); ++it) {
			vector<int> simplex_ids;
			for (auto it2 = (*it).S.begin(); it2 != (*it).S.end(); ++it2)
				simplex_ids.push_back(*it2);

			if (simplex_ids[0] < simplex_ids[1]) {
				D1(simplex_ids[0], (*it).id) = -1;
				D1(simplex_ids[1], (*it).id) = 1;
			}

			else {
				D1(simplex_ids[0], (*it).id) = 1;
				D1(simplex_ids[1], (*it).id) = -1;
			}
		}

		for (auto it = C[2].begin(); it != C[2].end(); ++it) {
			int highest_id = -1e6;
			int lowest_id = 1e6;
			for (auto it2 = (*it).S.begin(); it2 != (*it).S.end(); ++it2) {
				if ((*it2) > highest_id) highest_id = (*it2);
				if ((*it2) < lowest_id) lowest_id = (*it2);
			}

			simplex<int> temp;

			temp.S = (*it).S;
			temp.S.erase(highest_id);
			D2((*C[1].find(temp)).id, (*it).id) = 1;

			temp.S.clear();
			temp.S = (*it).S;
			temp.S.erase(lowest_id);
			D2((*C[1].find(temp)).id, (*it).id) = 1;

			temp.S.clear();
			temp.S.insert(lowest_id);
			temp.S.insert(highest_id);
			D2((*C[1].find(temp)).id, (*it).id) = -1;
		}

		arma::sp_mat L = (D1.t() * D1) + (D2 * D2.t());
//		arma::Mat<double> L = (D1.t() * D1) + (D2 * D2.t());
		arma::Mat<double> I = arma::eye(C[1].size(),C[1].size());
		arma::Mat<double> x = arma::zeros(C[1].size(),1);
		arma::Mat<double> x0 = arma::zeros(C[1].size(),1);
		arma::Mat<double> z = arma::zeros(C[2].size(),1);
		arma::Mat<double> f = arma::zeros(C[1].size(),1);
		arma::Mat<double> Lx;

		Mat homology_image_initial =
			Mat(ref_image.rows, ref_image.cols, CV_8UC3, Scalar (255, 255, 255));
		Mat homology_image_inter =
			Mat(ref_image.rows, ref_image.cols, CV_8UC3, Scalar (255, 255, 255));
		Mat homology_image_final =
			Mat(ref_image.rows, ref_image.cols, CV_8UC3, Scalar (255, 255, 255));

		vector<int> Avgs;
		for (int it = 0; it < C[1].size(); ++it) {
			int sum(0);
			for(auto p = id_to_simplex[it].S.begin(); p != id_to_simplex[it].S.end(); ++p)
				sum += landmarks[*p].obs_count;
			Avgs.push_back(sum);
		}

		int Avgs_Max = *max_element(Avgs.begin(), Avgs.end());
		for (int it = 0; it < Avgs.size(); ++it)
			x0(it,0) = exp(-double(Avgs[it])/Avgs_Max);

		if (x00.n_elem != 0) {
			for (int i = 0; i < x00.n_elem; ++i) {
				if (x00(i,0) == 0)
					x(i,0) = 1;

				else
					x(i,0) = 0.5;
			}

			for (int i2 = x00.n_elem; i2 < C[1].size(); ++i2)
				x(i2,0) = 1;
		}

		else
			x = x.col(0) + 1;

/*		draw_vector (x, homology_image_initial);
		imshow("homology_image_initial", homology_image_initial);
		waitKey (2e3);
*/
		double dt1 = 5e-2;

		do {
			Lx = L * x;
			x = x - Lx * dt1;
//			cout << "arma::norm(Lx): " << arma::norm(Lx) << endl;
		} while (arma::norm(Lx) > 1e-7);

		draw_vector (x, homology_image_inter);
		imshow("homology_image_inter", homology_image_inter);
		waitKey (2e3);

/*		do {
			f = ((D2.t() * D2) * z) + (2 * x.t() * D2).t();
			z = z - f * dt;
		} while (arma::norm(f) > 1e-4);
*/

//		x = (I  -  pinv(L) * L) * x;

		x = x / arma::norm(x,1);

		bool converged = false;
		double NormAvg(0), K_NormAvg(0);
		double NormAvgNext(0), K_NormAvgNext(0);
		int counter = 0;

		double dt2 = 5e-5;

		while (!converged) {
			f = D2*z + x;
			for (int i = 0; i < C[1].size(); ++i) {
				if (f(i,0) > 1e-7)
					f(i,0) = 1;

				else if (f(i,0) < -1e-7)
					f(i,0) = -1;

				else 
					f(i,0) = 0;
			}

			arma::Mat<double> sgn = D2.t() * f;

			z = z - sgn * dt2;

			NormAvg += arma::norm(sgn);
			counter++;
			if (counter == 100) {
				NormAvg = NormAvg/100;
				if (abs(NormAvg - NormAvgNext) < 1e-2) {
					converged = true;
				}

				else {
					NormAvgNext = NormAvg;
					NormAvg = 0;
					counter = 0;
				}
			}

//			cout << "arma::norm(x + D2*z,1): " << arma::norm(x + D2*z,1) << endl;
		}

		arma::Mat<double> K = x + D2 * z;
		H.clear();

		draw_vector (K, homology_image_final);
		imshow("homology_image_final", homology_image_final);
		waitKey (2e3);

		homology_image = nice_image.clone();

		for (int it = 0; it < C[1].size(); ++it) {
			double n = K(it,0);
			H.push_back (edges(abs(n), it));
		}

		sort(H.begin(), H.end(), edges_Comparator);

		vector<double> H_SD;
		for (int i = 0; i < H.size(); ++i)
			H_SD.push_back(H[i].homology_value);

		jump_pos = locate_Jump(H_SD);

		cout << "jump_pos: " << jump_pos << endl;

		x00 = arma::zeros(C[1].size(),1);
		for (int it = H.size()-1; it >= jump_pos; --it) {
			vector<int> C0_ids;
			x00(H[it].edge_id,0) = 1;
			for(auto p = id_to_simplex[H[it].edge_id].S.begin();
			p != id_to_simplex[H[it].edge_id].S.end(); ++p) {
				cc_ids.insert(*p);
				C0_ids.push_back(*p);
			}

			line (homology_image, landmarks[C0_ids[0]].coord,
				landmarks[C0_ids[1]].coord, Scalar (0, 0, 128), 1, 1);
		}

		identify_clusters();

		for (int it = 0; it < clusters.size(); ++it) {
			int OBSlnd(0), lnd(0);
			for (int it2 = 0; it2 < clusters[it].size(); ++it2) {
				if (clusters[it][it2].NearOBS)
					OBSlnd++;

				else
					lnd++;
			}

			if (lnd+OBSlnd > 3 && OBSlnd/(lnd+OBSlnd) >= 0.3) {
				for (int it2 = 0; it2 < clusters[it].size(); ++it2)
					cc_ids.erase(clusters[it][it2].id);

				clusters_ids.erase(clusters_ids.begin() + it);
			}
		}

		Mat FixLNDS(ref_image.rows, ref_image.cols, CV_8UC3, Scalar (255, 255, 255));
		for (int it = 0; it < clusters_ids.size(); ++it) {
			for (auto k = clusters_ids[it].begin(); k != clusters_ids[it].end(); ++k) {
				Point n = landmarks[*k].coord;
				for (int i = -1; i <= 1; ++i)
					for (int i2 = -1; i2 <= 1; ++i2)
						if (n.y+i >= 0 && n.y+i < FixLNDS.rows
						&& n.x+i2 >= 0 && n.x+i2 < FixLNDS.cols)
							FixLNDS.at<Vec3b> ((n.y)+i, (n.x)+i2) = Vec3b (10*it, 50*it, 100*it);
			}
		}

		imshow ("FixLNDS", FixLNDS);
		waitKey(1e2);

		cout << landmarks.size() << ", " << cc_ids.size() << ", " << H.size() << endl;

		char homology_image_name [50];
		sprintf (homology_image_name,
			"Results/Homology Image/homology_image_%d.png", homology_step);
		imwrite(homology_image_name, homology_image);
		imshow("homology_image", homology_image);
		waitKey (2e3);
	}

	void draw_simplicial_complex (Mat &image1, vector<Landmark> landmarks) {

		if ((C.size() > 1) && (C[1].size() != 0)) {
			for (auto it = C[1].begin(); it != C[1].end(); ++it) {
				for (auto it2 = (*it).S.begin(); it2 != (*it).S.end(); ++it2) {
					for (auto it3 = (*it).S.begin(); it3 != (*it).S.end(); ++it3) 
						line (image1, landmarks[*it2].coord, landmarks[*it3].coord,
						Scalar (0, 0, 0), 1, 1);
				}
			}
		}

		if ((C.size() > 2) && (C[2].size() != 0)) {
			for (auto it = C[2].begin(); it != C[2].end(); ++it) {
				vector<Point> elm;
				for (auto it2 = (*it).S.begin(); it2 != (*it).S.end(); ++it2)
					elm.push_back (landmarks[*it2].coord);
				const Point* ppt[1] = { &elm[0] };
				int npt[] = { 3 };
				Mat image2 = image1.clone();
				fillPoly (image1, ppt, npt, 1, Scalar (rand()%256, rand()%256 , rand()%256), 8);
				image1 = 0.6*image2 + 0.4*image1;
			}
		}
	}

	void display (int n) {

		if (n == 0) {
			for (auto it = C0.S.begin(); it != C0.S.end(); ++it)
				cout << (*it) << endl;
		}

		else {

			if (C.size() <= n) C.resize(n-1);
			if (n < C.size()) {
				for (auto it = C[n].begin(); it != C[n].end(); ++it) {
					for (auto it2 = (*it).S.begin(); it2 != (*it).S.end(); ++it2)
						cout << (*it2) << " ";
				cout << endl;
				}
			}

			else cout << "Simplex does not exist!" << endl;
		}
	}
};

simplicial_complex<int> sc;

class searchProblem: public AStar::Algorithm<Landmark,double> {
public:

	int rbt_number;
	vector<vector<Landmark>>* start_landmarks;
	vector<unordered_set<int>>* obs_count_list;
	vector<int> least_observed_count;
	vector<Landmark> least_observed_landmark;

	searchProblem (int rbt) {
		rbt_number = rbt;
		least_observed_count.resize(rbt);
		least_observed_landmark.resize(rbt);
		fill (least_observed_count.begin(), least_observed_count.end(), 1e6);
	}

	void getSuccessors (Landmark &n, std::vector<Landmark>* s, std::vector<double>* c) {
		for (auto i = n.neighbors.begin(); i != n.neighbors.end(); ++i) {
			if (n.neighbors_occlusion_count[*i] == 0) {
				(**i).rbt_id = n.rbt_id;
				s->push_back(**i);
				c->push_back((double)1.00);
			}
		}
	}

	std::vector<Landmark> getStartNodes (void) {
		std::vector<Landmark> startNodes;
		for (int r2 = 0; r2 < rbt_number; ++r2) {
			for (int a = 0; a < (*start_landmarks)[r2].size(); ++a) {
				startNodes.push_back ((*start_landmarks)[r2][a]);
			}
		}

		return (startNodes);
	}

	void nodeEvent (Landmark &n, unsigned int e) {
		if (e & EXPANDED) {
			landmarks[n.id].g_score = n.G;
		}
	}

	bool stopSearch (Landmark &n) {
		if (n.obs_count < least_observed_count[n.rbt_id]) {
			least_observed_landmark[n.rbt_id] = n;
			least_observed_count[n.rbt_id] = n.obs_count;
		}

		return false;
	}
};

class connected_components: public AStar::Algorithm<Landmark,double> {
public:

	Landmark* start_landmark;
	vector<Landmark> branch_landmarks;
	unordered_set<int> branch_landmarks_ids;

	void getSuccessors (Landmark &n, std::vector<Landmark>* s, std::vector<double>* c) {
		for (auto i = n.neighbors.begin(); i != n.neighbors.end(); ++i) {
			if (cc_ids.find((**i).id) != cc_ids.end()) {
				for (int it = H.size()-1; it >= jump_pos; --it) {
					vector<int> C0_ids;
					for(auto p = sc.id_to_simplex[H[it].edge_id].S.begin();
					p != sc.id_to_simplex[H[it].edge_id].S.end(); ++p)
						C0_ids.push_back(*p);
					if ((C0_ids[0] == n.id && C0_ids[1] == (**i).id)
					|| (C0_ids[1] == n.id && C0_ids[0] == (**i).id)) {
						s->push_back(**i);
						c->push_back((double)1.00);
					}
				}
			}
		}
	}

	std::vector<Landmark> getStartNodes (void) {
		std::vector<Landmark> startNodes;
		startNodes.push_back (*start_landmark);
		return (startNodes);
	}

	void nodeEvent (Landmark &n, unsigned int e) {
		if (e & EXPANDED) {
			branch_landmarks.push_back(n);
			branch_landmarks_ids.insert(n.id);
		}
	}
};

class HIW_search: public AStar::Algorithm<Landmark,double> {
public:

	vector<Landmark>* start_landmarks;
	Landmark nearest_landmark;

	void getSuccessors (Landmark &n, std::vector<Landmark>* s, std::vector<double>* c) {
		for (auto i = n.neighbors.begin(); i != n.neighbors.end(); ++i) {
			if (n.neighbors_occlusion_count[*i] == 0) {
				s->push_back(**i);
				c->push_back((double)1.00);
			}
		}
	}

	std::vector<Landmark> getStartNodes (void) {
		std::vector<Landmark> startNodes;
		for (int a = 0; a < (*start_landmarks).size(); ++a) {
			startNodes.push_back ((*start_landmarks)[a]);
		}

		return (startNodes);
	}

	bool stopSearch (Landmark &n) {
		if (cc_ids.find(n.id) != cc_ids.end()) {
			nearest_landmark = n;
			return true;
		}
	}
};

class Navigation: public AStar::Algorithm<Landmark,double> {
public:

	vector<Landmark>* start_landmarks;
	Landmark goal_landmark;

	Navigation (Landmark goal) {
		goal_landmark = goal;
	}

	void getSuccessors (Landmark &n, std::vector<Landmark>* s, std::vector<double>* c) {
		for (auto i = n.neighbors.begin(); i != n.neighbors.end(); ++i) {
			if (n.neighbors_occlusion_count[*i] == 0) {
				s->push_back(**i);
				c->push_back((double)1.00);
			}
		}
	}

	std::vector<Landmark> getStartNodes (void) {
		std::vector<Landmark> startNodes;
		for (int a = 0; a < (*start_landmarks).size(); ++a) {
			startNodes.push_back ((*start_landmarks)[a]);
		}

		return (startNodes);
	}

	bool stopSearch (Landmark &n) {
		return (n == goal_landmark);
	}
};

class robot_to_clusters: public AStar::Algorithm<Landmark,double> {
public:

	vector<Landmark>* start_landmarks;
	vector<unordered_set<int>> clusters_ids_copy;
	int robot_id;

	robot_to_clusters (int rbt_number) {
		robot_id = rbt_number;
		clusters_ids_copy = clusters_ids;
		observations[rbt_number].gscore_to_clusters.assign(clusters_ids_copy.size(), 1e6);
		observations[rbt_number].nearest_landmark_ids.clear();
		observations[rbt_number].nearest_landmark_ids.resize(clusters_ids_copy.size());
	}

	void getSuccessors (Landmark &n, std::vector<Landmark>* s, std::vector<double>* c) {
		for (auto i = n.neighbors.begin(); i != n.neighbors.end(); ++i) {
			s->push_back(**i);
			c->push_back((double)1.00);
		}
	}

	std::vector<Landmark> getStartNodes (void) {
		std::vector<Landmark> startNodes;
		for (int a = 0; a < (*start_landmarks).size(); ++a) {
			startNodes.push_back ((*start_landmarks)[a]);
		}

		return (startNodes);
	}

	void nodeEvent (Landmark &n, unsigned int e) {
		if (e & EXPANDED) {
			for (int it = 0; it < clusters_ids_copy.size(); ++it) {
				if (clusters_ids_copy[it].find(n.id) != clusters_ids_copy[it].end()
				&& observations[robot_id].gscore_to_clusters[it] == 1e6) {
					observations[robot_id].gscore_to_clusters[it] = n.G;
					observations[robot_id].nearest_landmark_ids[it] = n.id;
				}
			}
		}
	}
};
