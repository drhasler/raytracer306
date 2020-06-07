
// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(
        const char* fname,
        const std::vector<ConvexHull> &polygons,
        const std::vector<vec2> &seeds
        ) {
	FILE* f = fopen(fname, "w+");
	const char* cols[] = {
		"cornflowerblue",
		"darkslategray",
		"darkolivegreen",
		"lightcoral",
		"lightseagreen",
		"mediumturquoise",
		"olive",
		"rosybrown",
        "palevioletred",
	};
	int nc = sizeof(cols)/sizeof(char*);
    int c_idx = 0;
	fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
	for (const auto& pol: polygons) {
		fprintf(f, "<polygon points = \""); 
		for (vec2 p: pol.hull) {
			fprintf(f, "%3.3f, %3.3f ", p.x * 1000, 1000 - p.y * 1000);
		}
		fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n",
                cols[c_idx++ % nc]);
	}
    for (vec2 v: seeds)
		fprintf(f, "<circle cx=\"%3.3f\" cy=\"%3.3f\" r=\"4\" fill=\"white\""
                "/>\n", v.x * 1000, 1000 - v.y * 1000);
	fprintf(f, "</svg>\n");
	fclose(f);
}

#if 0
// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon {  
public:
	std::vector<Vector> vertices;
};	

// Adds one frame of an animated svg file. frameid is the frame number (between 0 and nbframes-1).
// polygons is a list of polygons, describing the current frame.
// The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
	void save_svg_animated(const std::vector<Polygon> &polygons, std::string filename, int frameid, int nbframes) {
		FILE* f;
		if (frameid == 0) {
			f = fopen(filename.c_str(), "w+");
			fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
			fprintf(f, "<g>\n");
		} else {
			f = fopen(filename.c_str(), "a+");
		}
		fprintf(f, "<g>\n");
		for (int i = 0; i < polygons.size(); i++) {
			fprintf(f, "<polygon points = \""); 
			for (int j = 0; j < polygons[i].vertices.size(); j++) {
				fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000-polygons[i].vertices[j][1] * 1000));
			}
			fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
		}
		fprintf(f, "<animate\n");
		fprintf(f, "	id = \"frame%u\"\n", frameid);
		fprintf(f, "	attributeName = \"display\"\n");
		fprintf(f, "	values = \"");
		for (int j = 0; j < nbframes; j++) {
			if (frameid == j) {
				fprintf(f, "inline");
			} else {
				fprintf(f, "none");
			}
			fprintf(f, ";");
		}
		fprintf(f, "none\"\n	keyTimes = \"");
		for (int j = 0; j < nbframes; j++) {
			fprintf(f, "%2.3f", j / (double)(nbframes));
			fprintf(f, ";");
		}
		fprintf(f, "1\"\n	dur = \"5s\"\n");
		fprintf(f, "	begin = \"0s\"\n");
		fprintf(f, "	repeatCount = \"indefinite\"/>\n");
		fprintf(f, "</g>\n");
		if (frameid == nbframes - 1) {
			fprintf(f, "</g>\n");
			fprintf(f, "</svg>\n");
		}
		fclose(f);
	}
#endif
