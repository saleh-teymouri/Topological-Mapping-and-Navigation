void RPGA (dpoint observations, int m,
		vector<dpoint> &Possible_Robot_Positions) { // Random Path Generation Algorithm

	int margin_cols = ref_image.cols/100;
	int margin_rows = ref_image.rows/100;
	int counter1(0);

	Start:
	double r = ((((double)rand())/((double)RAND_MAX))*(max_r-min_r)) + min_r;
//	max_d = PI * r;
	int intervals = max_d;
	double d = ((((double)rand())/((double)RAND_MAX))*(max_d-min_d)) + min_d;

	if (counter1 > 20)
		observations.angle = (((double)rand())/((double)RAND_MAX))*2*PI-PI;

	vector<Point> circle_pts;
	double alpha = observations.angle - (PI/2);
	double center_x = round(observations.coord.x + m*r*cos(alpha));
	double center_y = round(observations.coord.y + m*r*sin(alpha));
	Point center = Point(center_x, center_y);

	if (m == -1) {	/** m = -1 : counter-clockwise **/

		double start_angle = observations.angle - (PI/2);
		double end_angle = start_angle + (((double)d)/((double)r));
		double step_theta = ((double)(PI))/((double)(d)*(180));

		for (double theta = start_angle; theta <= end_angle; theta += step_theta)
			circle_pts.push_back(Point(round(center.x + (r*cos(theta))),
				round(center.y + (r*sin(theta)))));

		for (int i = 0; i < circle_pts.size(); ++i) {
			if (circle_pts[i].x >= ref_image.cols || circle_pts[i].x < 0
			|| circle_pts[i].y >= ref_image.rows || circle_pts[i].y < 0
			|| ref_image.at<Vec3b> (circle_pts[i].y, circle_pts[i].x) == obs_clr) {
				counter1++;
				goto Start;
			}
		}

		int it = circle_pts.size()-1;
		double new_x = circle_pts[it].x;
		double new_y = circle_pts[it].y;
		double new_direction = observations.angle + (((double)d)/((double)r));

		while (new_direction > PI)
			new_direction -= 2*PI;
		while (new_direction < -PI)
			new_direction += 2*PI;

		if (new_x < (ref_image.cols - margin_cols) && new_x >= (0 + margin_cols)
		&& new_y < (ref_image.rows - margin_rows) && new_y >= (0 + margin_rows)) {
			counter1 = 0;
			for (int h = 0; h < circle_pts.size(); h += (circle_pts.size()/intervals)) {
				double gamma = atan2 ((circle_pts[h].y - center.y), (circle_pts[h].x - center.x));
				Possible_Robot_Positions.push_back (dpoint (
					circle_pts[h].x, circle_pts[h].y, gamma+PI/2));
			}
		}

		else {
			counter1++;
			goto Start;
		}
	}

	if (m == 1) {	/** m = 1 : clockwise **/

		double start_angle = observations.angle + (PI/2);
		double end_angle = start_angle - (((double)d)/((double)r));
		double step_theta = ((double)(PI))/((double)(d)*(180));

		for (double theta = start_angle; theta >= end_angle; theta -= step_theta)
			circle_pts.push_back(Point(round(center.x + (r*cos(theta))),
				round(center.y + (r*sin(theta)))));

		for (int i = 0; i < circle_pts.size(); ++i) {
			if (circle_pts[i].x >= ref_image.cols || circle_pts[i].x < 0
			|| circle_pts[i].y >= ref_image.rows || circle_pts[i].y < 0
			|| ref_image.at<Vec3b> (circle_pts[i].y, circle_pts[i].x) == obs_clr) {
				counter1++;
				goto Start;
			}
		}

		int it = circle_pts.size()-1;
		double new_x = circle_pts[it].x;
		double new_y = circle_pts[it].y;
		double new_direction = observations.angle - (((double)d)/((double)r));

		while (new_direction > PI)
			new_direction -= 2*PI;
		while (new_direction < -PI)
			new_direction += 2*PI;

		if (new_x < (ref_image.cols - margin_cols) && new_x >= (0 + margin_cols)
		&& new_y < (ref_image.rows - margin_rows) && new_y >= (0 + margin_rows)) {
			counter1 = 0;
			for (int h = 0; h < circle_pts.size(); h += (circle_pts.size()/intervals)) {
				double gamma = atan2 ((circle_pts[h].y - center.y), (circle_pts[h].x - center.x));
				Possible_Robot_Positions.push_back (dpoint (
					circle_pts[h].x, circle_pts[h].y, gamma-PI/2));
			}
		}

		else {
			counter1++;
			goto Start;
		}
	}
}

void Random_Walk (int id) {

	vector<dpoint> Possible_Robot_Positions;
	double d = observations[id].sensor_radius;
	double m = pow (-1, rand()%2);

	RPGA (observations[id], m, Possible_Robot_Positions);

	while (Possible_Robot_Positions.size() != 0) {
		observations_mtx.lock();
		observations[id].coord = Possible_Robot_Positions[0].coord;
		observations[id].angle = Possible_Robot_Positions[0].angle;
		observations_mtx.unlock();
		Possible_Robot_Positions.erase (Possible_Robot_Positions.begin());
		exploration (id);

		if (Show_image) {
			disp (sc_image, observations[id].coord, Vec3b(128*id/2, 0, 255*id/2), 1);
//			Mat rbt_image = sc_image.clone();
			draw_robot (rbt_image, observations[id]);
//			imshow("img0",rbt_image);
//			waitKey(animation_wait);
		}
	}
}

void Navigate (int id, vector<Landmark> path_landmarks, bool & success) {

	double d = observations[id].sensor_radius;
	double angle_range = observations[id].sensor_angle;
	int m(0), counter1(0), counter2(0), counter3(0), counter4(0);

	Start:
	while (path_landmarks.size() != 0) {

		if (Show_image) {
			img_mtx.lock();
			disp (sc_image, observations[id].coord, Vec3b(128*id/2, 0, 255*id/2), 1);
			draw_path (PthHmlgy_image, id, path_landmarks, Scalar(128*id/2, 0, 255*id/2));
//			Mat rbt_image = sc_image.clone();
			for (int i = 0; i != path_landmarks.size(); ++i)
				disp(rbt_image, path_landmarks[i].coord, Vec3b(255, 0, 0), 2);
			draw_robot (rbt_image, observations[id]);
//			imshow("img0",rbt_image);
//			waitKey(animation_wait);
			img_mtx.unlock();
		}

		bool no_observation = true;

		for (int n1 = path_landmarks.size()-1; n1 >= 0; --n1) {

			path_landmarks[n1].distances = 
				(dist2(observations[id].coord, path_landmarks[n1].coord, d));
		
			if (path_landmarks[n1].distances < d
			&& line_check2(observations[id].coord,
				path_landmarks[n1].coord, ref_image, obs_clr)
			&& angle_check(observations[id], path_landmarks[n1].coord, angle_range)) {
				vector<dpoint> Possible_Robot_Positions;
				left_or_right (observations[id], path_landmarks[n1], m);
				RPGA (observations[id], m, Possible_Robot_Positions);

				while (Possible_Robot_Positions.size() != 0) {
					observations_mtx.lock();
					observations[id].coord = Possible_Robot_Positions[0].coord;
					observations[id].angle = Possible_Robot_Positions[0].angle;
					observations_mtx.unlock();
					Possible_Robot_Positions.erase (Possible_Robot_Positions.begin());
					left_or_right (observations[id], path_landmarks[n1], m);
					exploration (id);

					bool not_visible = true;

					for (int m1 = path_landmarks.size()-1; m1 >= 0; --m1) {
						path_landmarks[m1].distances = 
						(dist2(observations[id].coord, path_landmarks[m1].coord, d));

						if (path_landmarks[m1].distances < d
						&& line_check2(observations[id].coord,
							path_landmarks[m1].coord, ref_image, obs_clr)
						&& angle_check(observations[id],
							path_landmarks[m1].coord, angle_range)) {

							if (path_landmarks[m1].coord != 
								path_landmarks[n1].coord) {

								if (m1 == 0)
									path_landmarks.erase (path_landmarks.begin());

								else
									path_landmarks.erase (path_landmarks.begin(),
										path_landmarks.begin()+m1);

								goto Start;
							}

							if (path_landmarks.size() == 1)
								path_landmarks.erase (path_landmarks.begin());

							not_visible = false;
							counter1 = 0;
							counter3 = 0;
						}
					}

					if (not_visible) {
						observations_mtx.lock();
						if (m == -1)
							observations[id].angle += (0.9*angle_range);
						if (m == 1)
							observations[id].angle -= (0.9*angle_range);
						while (observations[id].angle > PI)
							observations[id].angle -= 2*PI;
						while (observations[id].angle < -PI)
							observations[id].angle += 2*PI;
						observations_mtx.unlock();
						path_landmarks.erase (path_landmarks.begin(),
							path_landmarks.begin()+n1);
						exploration (id);
						counter1++;
						if (counter1 > ((2*PI)/angle_range)) {
							if (counter3 < 1) {
								Random_Walk (id);
								counter3++;
							}

							else {
								success = false;
								if (path_landmarks.size() > 1) {
									landmarks[path_landmarks[0].id].neighbors_occlusion_count[
										&landmarks[path_landmarks[1].id]]++;
									landmarks[path_landmarks[1].id].neighbors_occlusion_count[
										&landmarks[path_landmarks[0].id]]++;
								}

								return;
							}
						}

						goto Start;
					}

					if (Show_image) {
						img_mtx.lock();
						disp (sc_image, observations[id].coord, Vec3b(128*id/2, 0, 255*id/2), 1);
						draw_path (PthHmlgy_image, id, path_landmarks,
							Scalar(128*id/2, 0, 255*id/2));
//						Mat rbt_image = sc_image.clone();
						for (int i = 0; i != path_landmarks.size(); ++i)
							disp(rbt_image, path_landmarks[i].coord, Vec3b(255, 0, 0), 2);
						draw_robot (rbt_image, observations[id]);
//						imshow("img0",rbt_image);
//						waitKey(animation_wait);
						img_mtx.unlock();
					}
				}

				if (path_landmarks.size() == 1)
					path_landmarks.erase (path_landmarks.begin());
				else
					path_landmarks.erase (
						path_landmarks.begin(), path_landmarks.begin()+n1);

				no_observation = false;
				counter2 = 0;
				counter4 = 0;
			}
		}

		if (no_observation) {
			observations_mtx.lock();
			if (m == -1)
				observations[id].angle += (0.9*angle_range);
			if (m == 1)
				observations[id].angle -= (0.9*angle_range);
			while (observations[id].angle > PI)
				observations[id].angle -= 2*PI;
			while (observations[id].angle < -PI)
				observations[id].angle += 2*PI;
			observations_mtx.unlock();
			exploration (id);
			counter2++;
			if (counter2 > ((2*PI)/angle_range)) {
				if (counter4 < 1) {
					Random_Walk (id);
					counter4++;
				}

				else {
					success = false;
					if (path_landmarks.size() > 1) {
						landmarks_mtx.lock();
						landmarks[path_landmarks[0].id].neighbors_occlusion_count[
							&landmarks[path_landmarks[1].id]]++;
						landmarks[path_landmarks[1].id].neighbors_occlusion_count[
							&landmarks[path_landmarks[0].id]]++;
						landmarks_mtx.unlock();
					}

					return;
				}
			}

			goto Start;
		}
	}
}

void Combined_RW_ISW (int id) {

	Start:

	vector<vector<Landmark>> all_observing_landmarks(rbt);

	exploration (id);

	observations_mtx.lock();
	for (int r1 = 0; r1 < rbt; ++r1)
		for (int m1 = 0; m1 < observations[r1].observing_landmarks.size(); ++m1)
			all_observing_landmarks[r1].push_back (observations[r1].observing_landmarks[m1]);

	searchProblem test_search(rbt);
	test_search.start_landmarks = &all_observing_landmarks;
	test_search.obs_count_list = &obs_list;
	test_search.search();
	observations_mtx.unlock();

//	cout << test_search.least_observed_count[id] << ", " <<
//		test_search.least_observed_landmark[id].obs_count << endl;

	if (test_search.least_observed_landmark[id].g_score < 5) {
		for (int k1 = 0; k1 < 5; ++k1)
			Random_Walk (id);
		ISW_counter = 0;
	}

/*	else {
		for (int k2 = 0; k2 < obs_list.size(); ++k2) {
			if (obs_list[k2].find(test_search.least_observed_landmark[id].id)
			!= obs_list[k2].end())
				cout << k2 << ", " << test_search.least_observed_count[id] << ", " <<
					test_search.least_observed_landmark[id].obs_count << endl;
		}
	}
*/
	vector<Landmark> path_landmarks;
	std::vector<Landmark*> path =
		test_search.reconstructPointerPath (test_search.least_observed_landmark[id]);
	for (int p2 = path.size()-1; p2 >= 0; --p2) {
		path_landmarks.push_back (Landmark (
			(*path[p2]).coord.x, (*path[p2]).coord.y, (*path[p2]).id));
	}

	bool success = true;
	Navigate (id, path_landmarks, success);

	int counter1(0);
	while (!success) {
		counter1++;
		Navigation find_path(test_search.least_observed_landmark[id]);
		find_path.start_landmarks = &observations[id].observing_landmarks;
		find_path.search();
		vector<Landmark> new_path_landmarks;
		std::vector<Landmark*> new_path =
			find_path.reconstructPointerPath (test_search.least_observed_landmark[id]);
		for (int p2 = new_path.size()-1; p2 >= 0; --p2) {
			new_path_landmarks.push_back (Landmark (
				(*new_path[p2]).coord.x, (*new_path[p2]).coord.y, (*new_path[p2]).id));
		}

		success = true;
		Navigate (id, new_path_landmarks, success);
		if (counter1 > 20)
			goto Start;
	}

	ISW_counter++;

	if (ISW_counter > factor*num_ISW) {
		for (int k3 = 0; k3 < factor*num_RW; ++k3)
			Random_Walk (id);
		ISW_counter = 0;
	}
}

void Informed_Systematic_Walk (int id) {

	Start:

	vector<vector<Landmark>> all_observing_landmarks(rbt);

	exploration (id);

	observations_mtx.lock();
	for (int r1 = 0; r1 < rbt; ++r1)
		for (int m1 = 0; m1 < observations[r1].observing_landmarks.size(); ++m1)
			all_observing_landmarks[r1].push_back (observations[r1].observing_landmarks[m1]);

	searchProblem test_search(rbt);
	test_search.start_landmarks = &all_observing_landmarks;
	test_search.obs_count_list = &obs_list;
	test_search.search();
	observations_mtx.unlock();

	vector<Landmark> path_landmarks;
	std::vector<Landmark*> path =
		test_search.reconstructPointerPath (test_search.least_observed_landmark[id]);
	for (int p2 = path.size()-1; p2 >= 0; --p2) {
		path_landmarks.push_back (Landmark (
			(*path[p2]).coord.x, (*path[p2]).coord.y, (*path[p2]).id));
	}

	bool success = true;
	Navigate (id, path_landmarks, success);

	int counter1(0);
	while (!success) {
		counter1++;
		Navigation find_path(test_search.least_observed_landmark[id]);
		find_path.start_landmarks = &observations[id].observing_landmarks;
		find_path.search();
		vector<Landmark> new_path_landmarks;
		std::vector<Landmark*> new_path =
			find_path.reconstructPointerPath (test_search.least_observed_landmark[id]);
		for (int p2 = new_path.size()-1; p2 >= 0; --p2) {
			new_path_landmarks.push_back (Landmark (
				(*new_path[p2]).coord.x, (*new_path[p2]).coord.y, (*new_path[p2]).id));
		}

		success = true;
		Navigate (id, new_path_landmarks, success);
		if (counter1 > 20)
			goto Start;
	}
}

void Homology_Informed_Walk (int id) {

	Start:

	if (cc_ids.size() > 0) {

		observations[id].Assigned = true;
		connected_components cc;
		bool success;
		do {
			observations_mtx.lock();
			HIW_search closest_cluster;
			closest_cluster.start_landmarks = &observations[id].observing_landmarks;
			closest_cluster.search();
			observations_mtx.unlock();

//			disp (clusters_img, closest_cluster.nearest_landmark.coord, Vec3b(255, 0, 0), 2);

			cc.start_landmark = &landmarks[closest_cluster.nearest_landmark.id];
			cc.search();

			vector<Landmark> path_landmarks;
			std::vector<Landmark*> path =
				closest_cluster.reconstructPointerPath (closest_cluster.nearest_landmark);
			for (int p2 = path.size()-1; p2 >= 0; --p2) {
				path_landmarks.push_back (Landmark (
					(*path[p2]).coord.x, (*path[p2]).coord.y, (*path[p2]).id));
			}

			success = true;
			Navigate (id, path_landmarks, success);

		} while (!success);

		for (int it2 = 0; it2 < cc.branch_landmarks.size(); ++it2) {
//			disp (clusters_img, cc.branch_landmarks[it2].coord, Vec3b(0, 0, 255), 1);
			cc_ids.erase(cc.branch_landmarks[it2].id);
		}

		while (cc.branch_landmarks.size() > 0) {
			int rand_it = ((((double)rand())/((double)RAND_MAX))*(cc.branch_landmarks.size()));
			int counter1(0);
			do {
				counter1++;
				Navigation HIW_exploration(cc.branch_landmarks[rand_it]);
				HIW_exploration.start_landmarks = &observations[id].observing_landmarks;
				HIW_exploration.search();

				vector<Landmark> HIW_path_landmarks;
				std::vector<Landmark*> HIW_path =
					HIW_exploration.reconstructPointerPath (cc.branch_landmarks[rand_it]);
				for (int p2 = HIW_path.size()-1; p2 >= 0; --p2) {
					HIW_path_landmarks.push_back (Landmark (
						(*HIW_path[p2]).coord.x, (*HIW_path[p2]).coord.y, (*HIW_path[p2]).id));
				}

				success = true;
				Navigate (id, HIW_path_landmarks, success);
				if (counter1 > 20)
					goto Start;

			} while (!success);

			cc.branch_landmarks.erase (cc.branch_landmarks.begin() + rand_it);
		}

		observations[id].Assigned = false;
	}

	else
		HIW_finished[id] = true;
}

void updated_Homology_Informed_Walk (int id) {

	Start:

	if (clusters_ids.size() > 0) {

		observations[id].Assigned = true;
		observations_mtx.lock();
		for (int r2 = 0; r2 < rbt; ++r2) {
			robot_to_clusters rbt2cluster (r2);
			rbt2cluster.start_landmarks = &observations[r2].observing_landmarks;
			rbt2cluster.search();
		}

		vector<vector<double>> costMatrix(rbt);
		for (int it1 = 0; it1 < rbt; ++it1)
			for (int it2 = 0; it2 < observations[it1].gscore_to_clusters.size(); ++it2)
				costMatrix[it1].push_back(observations[it1].gscore_to_clusters[it2]);

		HungarianAlgorithm HungAlgo;
		vector<int> assignment;

		double cost = HungAlgo.Solve(costMatrix, assignment);
		int AsgnClstr = assignment[id];

//		draw_homology_image(assignment);

		if (AsgnClstr < 0 || AsgnClstr > clusters_ids.size()) {
			HIW_finished[id] = true;
			observations[id].Assigned = false;
			observations_mtx.unlock();
			return;
		}

		unordered_set<int> Assigned_Cluster_ids = clusters_ids[AsgnClstr];
		vector<Landmark> Assigned_Cluster = clusters[AsgnClstr];

		int Nearest_Landmark_id = observations[id].nearest_landmark_ids[AsgnClstr];
		Landmark Nearest_Landmark = landmarks[Nearest_Landmark_id];

/*		if (Nearest_Landmark.coord.x == 0 && Nearest_Landmark.coord.y == 0) {
			cout << "ERROR: " << id << ", AsgnClstr: " << AsgnClstr << ", size: " <<
			observations[id].nearest_landmark_ids.size() << ", clusters_ids_size: " <<
			clusters_ids.size() << ", nrstlnd_id: " << Nearest_Landmark_id << endl;
		}
*/
		clusters_ids.erase(clusters_ids.begin() + AsgnClstr);
		clusters.erase(clusters.begin() + AsgnClstr);
		cc_ids.erase(Nearest_Landmark_id);
		for (int i = 0; i < Assigned_Cluster.size(); ++i)
			cc_ids.erase(Assigned_Cluster[i].id);

/*		for (unsigned int x = 0; x < costMatrix.size(); x++)
			cout << x << "," << assignment[x] << "\t";

		cout << "\ncost: " << cost << endl;

		for (int it = 0; it < rbt; ++it) {
			cout << it << "\t";
			for (int it2 = 0; it2 < observations[it].gscore_to_clusters.size(); ++it2)
				cout << observations[it].gscore_to_clusters[it2] << ", ";
			cout << endl;
		}
*/

		observations_mtx.unlock();

		bool success;
		do {
			observations_mtx.lock();
			Navigation HIW_exploration(Nearest_Landmark);
			HIW_exploration.start_landmarks = &observations[id].observing_landmarks;
			HIW_exploration.search();
			observations_mtx.unlock();

			vector<Landmark> HIW_path_landmarks;
			std::vector<Landmark*> HIW_path =
				HIW_exploration.reconstructPointerPath(Nearest_Landmark);
			for (int p2 = HIW_path.size()-1; p2 >= 0; --p2) {
				HIW_path_landmarks.push_back (Landmark (
					(*HIW_path[p2]).coord.x, (*HIW_path[p2]).coord.y, (*HIW_path[p2]).id));
			}

			success = true;
			if (HIW_path_landmarks.size() > 1)
				Navigate (id, HIW_path_landmarks, success);

/*			else {
				cout << "\nrbt_p "<<id<<": "<< observations[id].coord <<
				", nrstlnd_p: " << Nearest_Landmark.coord << ", id: " <<
				Nearest_Landmark_id << endl;
			}
*/
		} while (!success);

		while (Assigned_Cluster.size() > 0) {
			int rand_it = ((((double)rand())/((double)RAND_MAX))*(Assigned_Cluster.size()));
			int counter1(0);
			do {
				counter1++;
				observations_mtx.lock();
				Navigation HIW_exploration(Assigned_Cluster[rand_it]);
				HIW_exploration.start_landmarks = &observations[id].observing_landmarks;
				HIW_exploration.search();
				observations_mtx.unlock();

				vector<Landmark> HIW_path_landmarks;
				std::vector<Landmark*> HIW_path =
					HIW_exploration.reconstructPointerPath (Assigned_Cluster[rand_it]);
				for (int p2 = HIW_path.size()-1; p2 >= 0; --p2) {
					HIW_path_landmarks.push_back (Landmark (
						(*HIW_path[p2]).coord.x, (*HIW_path[p2]).coord.y, (*HIW_path[p2]).id));
				}

				success = true;
				Navigate (id, HIW_path_landmarks, success);
				if (counter1 > 20)
					goto Start;

			} while (!success);

			Assigned_Cluster.erase(Assigned_Cluster.begin() + rand_it);

		}

		observations[id].Assigned = false;
	}

	else
		HIW_finished[id] = true;
}

void move_robot(int id) {

	while (sc.C0.S.size() < landmarks.size()/10)
		Random_Walk (observations[id].rbt_id);

	while ((double(C2_size)/sc_number) <= first_percent_complete) //rt > 0.05
		Combined_RW_ISW (observations[id].rbt_id);

/*	while ((double(C2_size)/sc_number) <= second_percent_complete)
		Informed_Systematic_Walk (observations[id].rbt_id);
*/
	RWISW_finished[id] = true;
}
