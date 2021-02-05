bool line_check (Point pt1, Point pt2, Mat image) {

	vector<int> color;
	vector<Point> pts;
	
	line_coords (pt1, pt2, pts);

	for (int i = 0; i != pts.size(); ++i)
		color.push_back (image.at<Vec3b> ((pts[i].y), pts[i].x)[0]);

	for (int m = 0; m != color.size(); ++m)
		if (color[m] != 255) return false;
}

bool line_check2 (Point pt1, Point pt2, Mat image, Vec3b c) {

	Point pts;
	double big(1);

	if (abs(pt2.x - pt1.x) > abs(pt2.y - pt1.y)) { big = abs(pt2.x - pt1.x); }
	else { big = abs(pt2.y - pt1.y); }

	for (double n = 0; n <= abs(pt2.x - pt1.x) || n <= abs(pt2.y - pt1.y); ++n) {
		pts = (Point(
		(pt1.x + round((double)(((double)n)*((double)(pt2.x - pt1.x)))/((double)big))),
		(pt1.y + round((double)(((double)n)*((double)(pt2.y - pt1.y)))/((double)big)))));

		if (pts.x >= image.cols || pts.x < 0 || pts.y >= image.rows || pts.y < 0
		|| image.at<Vec3b> (pts.y, pts.x) == c)
			return false;
	}
	return true;
}

bool angle_check (dpoint pt1, Point pt2, double angle_range) {

	double angle(0);

	angle = atan2 ((pt2.y - pt1.coord.y), (pt2.x - pt1.coord.x));

	if (abs(angle - pt1.angle) < (angle_range/2)
	|| abs(angle - pt1.angle + 2*PI) < (angle_range/2)
	|| abs(angle - pt1.angle - 2*PI) < (angle_range/2))
		return true;
	else 
		return false;
}

bool path_check (Point center, double r, Mat image) {

	Point pts;
	int color;

	for (double theta = 0; theta < 2*PI; theta += 0.01) {
		pts = (Point((center.x + (r*cos(theta))), (center.y + (r*sin(theta)))));
		if (pts.x > image.cols || pts.x < 0 || pts.y > image.rows || pts.y < 0)
			return false;
		else {
			color = (image.at<Vec3b> ((pts.y), pts.x)[0]);
			if (color == 0)
				return false;
		}
	}
	return true;
}

bool boundary_check (Point pts, Mat image, Vec3b c) {

	for (int r = 0; r <= 1; ++r) {
		for (int r2 = 0; r2 <= 1; ++r2) {
			if (pts.x+r2 >= image.cols || pts.x+r2 < 0 || pts.y+r >= image.rows || pts.y+r < 0
			|| image.at<Vec3b> (pts.y+r, pts.x+r2) == c)
				return false;
		}
	}
	return true;
}

bool validation_check (Point pt1, Point pt2, Mat image, Vec3b c) {

	Point pts;
	double big(1);

	if (abs(pt2.x - pt1.x) > abs(pt2.y - pt1.y)) { big = abs(pt2.x - pt1.x); }
	else { big = abs(pt2.y - pt1.y); }

	for (double n = 0; n <= abs(pt2.x - pt1.x) || n <= abs(pt2.y - pt1.y); ++n) {
		pts = (Point(
		(pt1.x + round((double)(((double)n)*((double)(pt2.x - pt1.x)))/((double)big))),
		(pt1.y + round((double)(((double)n)*((double)(pt2.y - pt1.y)))/((double)big)))));

		for(int r = 0; r <= 1; ++r) {
			for(int r2 = 0; r2 <= 1; ++r2) {
				if (pts.x+r2 >= image.cols || pts.x+r2 < 0 || pts.y+r >= image.rows || pts.y+r < 0
				|| image.at<Vec3b> (pts.y+r, pts.x+r2) == c)
					return false;
			}
		}
	}
	return true;
}

void exploration (int id) {

	expl_count++;

	simplex<int> landmarks_id;
	vector<Landmark> sorted_landmarks;
	double d = observations[id].sensor_radius;
	double angle_range = observations[id].sensor_angle;

	start:
	for (int i = 0; i != landmarks.size(); ++i) {
		landmarks_mtx.lock();
		landmarks[i].distances = (dist2(observations[id].coord, landmarks[i].coord, d));
		if (landmarks[i].distances < d
		&& line_check2(observations[id].coord, landmarks[i].coord, ref_image, obs_clr)
		&& angle_check(observations[id], landmarks[i].coord, angle_range)) {
			landmarks[i].rbt_id = id;
			sorted_landmarks.push_back (landmarks[i]);
			for (int n = -3; n <= 3; ++n) {
				bool BRK = false;
				for (int m = -3; m <= 3; ++m) {
					Point P1 = Point(landmarks[i].coord.x+n, landmarks[i].coord.y+m);
					double Dstnc = dist2(observations[id].coord, P1, d);
					if (Dstnc < d && P1.y >= 0 && P1.y < ref_image.rows
					&& P1.x >= 0 && P1.x < ref_image.cols
					&& ref_image.at<Vec3b> (P1.y, P1.x) == obs_clr) {
						landmarks[i].NearOBS = true;
						BRK = true;
						break;
					}
				}

				if (BRK)
					break;
			}
		}

		landmarks_mtx.unlock();
	}

	observations_mtx.lock();
	if (sorted_landmarks.empty()) {
		observations[id].angle = (((double)rand())/((double)RAND_MAX))*2*PI-PI;
		observations_mtx.unlock();
		goto start;
	}

	if (!observations[id].observing_landmarks.empty()) {
		observations[id].observing_landmarks.erase
			(observations[id].observing_landmarks.begin(),
			observations[id].observing_landmarks.end());
	}

	sort (sorted_landmarks.begin(), sorted_landmarks.end(), Landmark_Comparator);
	observations[id].observing_landmarks = sorted_landmarks;
	observations_mtx.unlock();

	for (int k1 = 0; k1 < sorted_landmarks.size(); ++k1) {
		for (int k2 = 0; k2 < sorted_landmarks.size(); ++k2) {
			landmarks_mtx.lock();
			landmarks[sorted_landmarks[k1].id].neighbors.insert
				(&landmarks[sorted_landmarks[k2].id]);
			landmarks_mtx.unlock();
		}
	}

	for (int i2 = 0; i2 < 7 && i2 < sorted_landmarks.size(); ++i2) {
		landmarks_id.S.insert (sorted_landmarks[i2].id);
		landmarks_mtx.lock();
		landmarks[sorted_landmarks[i2].id].obs_count++;
		if (landmarks[sorted_landmarks[i2].id].obs_count + 1 > obs_list.size())
			obs_list.resize(landmarks[sorted_landmarks[i2].id].obs_count + 1);
		obs_list[landmarks[sorted_landmarks[i2].id].obs_count - 1].erase(sorted_landmarks[i2].id);
		obs_list[landmarks[sorted_landmarks[i2].id].obs_count].insert(sorted_landmarks[i2].id);
		landmarks_mtx.unlock();
	}

	sc_mtx.lock();
	sc.c2_pointers.clear();
	sc.create_observation (landmarks_id);
	sc.redundant_simplex (landmarks_id);
	if (live_simplical_comlex_image)
		draw_simplex (simplex_image, sorted_landmarks);
	sc_mtx.unlock();
}

void RS_exploration (int id) {

	expl_count++;

	simplex<int> landmarks_id;
	vector<Landmark> sorted_landmarks;
	double d = observations[id].sensor_radius;
	double angle_range = observations[id].sensor_angle;

	for (int i = 0; i != landmarks.size(); ++i) {
		landmarks[i].distances = (dist(observations[id].coord, landmarks[i].coord));
		if (landmarks[i].distances < d
		&& line_check2(observations[id].coord, landmarks[i].coord, ref_image, obs_clr)
		&& angle_check(observations[id], landmarks[i].coord, angle_range))
			sorted_landmarks.push_back (landmarks[i]);
	}
	sort (sorted_landmarks.begin(), sorted_landmarks.end(), Landmark_Comparator);
	for (int i2 = 0; i2 < 7 && i2 < sorted_landmarks.size(); ++i2)
		landmarks_id.S.insert (sorted_landmarks[i2].id);

	sc.c2_pointers.clear();
	sc.create_observation (landmarks_id);
	sc.redundant_simplex (landmarks_id);
}

void Raster_scan () {

	for (int j = 0; j < ref_image.rows; ++j) {
		for (int i = 0; i < ref_image.cols; ++i) {
			for (double theta = -PI; theta < PI; theta += PI/30) {
				cout<<"expl_count: "<<expl_count<<"\t\tsc_count: "<<sc_count
				<<"\t\t\t\tc[2].size: "<<C2_size<<"\r"<<flush;
				if (ref_image.at<Vec3b> (j, i) != obs_clr) {
					observations[0].coord = Point (i, j);
					observations[0].angle = theta;
					RS_exploration (observations[0].rbt_id);
				}
			}
		}
	}
}
