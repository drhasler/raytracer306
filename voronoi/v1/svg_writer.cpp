const char* nice_colors[] = {
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

void save_svg(const char* fname, const Diagram& diagram) {
	FILE* f = fopen(fname, "w+");
	int nc = sizeof(nice_colors)/sizeof(char*);
    int cnt = 0;
	fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
	for (const auto& cell: diagram) {
		fprintf(f, "<polygon points=\""); 
		for (vec2 p: cell.pts) {
            p += cell.kern;
			fprintf(f, "%3.3f, %3.3f ", p.x * 1000, 1000 - p.y * 1000);
		}
		fprintf(f, "\"\nfill=\"%s\" stroke=\"black\"/>\n",
                nice_colors[cnt++ % nc]);
	}
	for (const auto& cell: diagram)
		fprintf(f, "<circle cx=\"%3.3f\" cy=\"%3.3f\" r=\"4\" fill=\"white\""
                "/>\n", cell.kern.x * 1000, 1000 - cell.kern.y * 1000);
	fprintf(f, "</svg>\n");
	fclose(f);
}
