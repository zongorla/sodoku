#include <cstdio>
#include <string>
#include <vector>
#include <vector>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <unistd.h>
typedef unsigned int uint;
class Solver
{
public:
	Solver();
	Solver(std::string init, char zerochar);
	Solver(const Solver* init);

	~Solver();

	void print(std::ostream &s);

	bool isSolved();
	bool isAllowed(char val, int x, int y);
	bool solveBackTrack();

	void set(char val, int x, int y);

private:
	int width;
	int fullsize;
	int block;
	char zerochar;
	char **data;

};

Solver::Solver()
{
	for (int y = 0; y < 9; ++y)
	{
		for (int x = 0; x < 9; ++x)
		{
			data[y][x] = 0;
		}
	}
}

Solver::Solver(std::string init, char zerochararg)
{
	fullsize = init.length();
	width = sqrt(fullsize);
	block = sqrt(width);
	zerochar = zerochararg;
	data = (char**)malloc(sizeof(char*) * width);
	for (int i = 0; i < width; i++) {
		data[i] = (char*)malloc(sizeof(char)* width);
	}
	for (int i = 0; i < fullsize; ++i)
	{
		int x = i % width;
		int y = i / width;
		data[y][x] = init[i] - zerochar;
	}
}
Solver::Solver(const Solver * init)
{
	fullsize = init->fullsize;
	width = sqrt(fullsize);
	block = sqrt(width);
	zerochar = init->zerochar;
	data = (char**)malloc(sizeof(char*) * width);
	for (int i = 0; i < width; i++) {
		data[i] = (char*)malloc(sizeof(char)* width);
	}
	for (int y = 0; y < init->width; ++y)
	{
		for (int x = 0; x < init->width; ++x)
		{
			data[y][x] = init->data[y][x];
		}
	}
}


Solver::~Solver()
{
}

void Solver::print(std::ostream & s)
{
	for (int y = 0; y < width; ++y)
	{
		for (int x = 0; x < width; ++x)
		{
			s << (char)(data[y][x] + zerochar) << " ";
		}
		s << std::endl;
	}
}

bool Solver::isSolved()
{
	// Minden cella ki van töltve a táblában?
	for (int y = 0; y < width; ++y)
	{
		for (int x = 0; x < width; ++x)
		{
			if (data[y][x] == 0) return false;
		}
	}
	return true;
}

bool Solver::isAllowed(char val, int x, int y)
{
	bool allowed = true;

	// Azonos sorban vagy oszlopban csak egy 'val' lehet
	for (int i = 0; i < width; ++i)
	{
		if (data[y][i] == val) allowed = false;
		if (data[i][x] == val) allowed = false;
	}

	// Az adott blockxblock-as cellában csak egy 'val' lehet
	int cellBaseX = block * (int)(x / block);
	int cellBaseY = block * (int)(y / block);
	for (int y = cellBaseY; y < cellBaseY + block; ++y)
	{
		for (int x = cellBaseX; x < cellBaseX + block; ++x)
		{
			if (data[y][x] == val) allowed = false;
		}
	}

	return allowed;
}

bool Solver::solveBackTrack()
{
	// Készen vagyunk?
	if (isSolved())
	{
		return true;
	}

	// Keressünk egy pozíciót, amely még nincs kitöltve
	for (int y = 0; y < width; ++y)
	{
		for (int x = 0; x < width; ++x)
		{
			// Nincs még kitöltve?
			if (data[y][x] == 0)
			{
				// Keressünk egy értéket, amely megfelel a szabályoknak
				for (int n = 1; n <= width; ++n)
				{
					// Beírható az adott pozícióba?
					if (isAllowed(n, x, y))
					{
						// Másoljuk le a táblát
						Solver tmpSolver(this);
						// Írjuk bele az új értéket
						tmpSolver.set(n, x, y);
						// Próbáljuk megoldani az új táblát
						if (tmpSolver.solveBackTrack())
						{
							// Megoldás
							*this = tmpSolver;
							return true;
						}
					}
				}
			}
			// Nem tudtunk értéket írni a cellába, így lépjünk vissza
			if (data[y][x] == 0) return false;
		}
	}

	return false;
}

void Solver::set(char val, int x, int y)
{
	data[y][x] = val;
}


using namespace std;
bool debugtimers = true;
class TimeInterval {


	clock_t _start;
	clock_t _end;
	bool started = false;
	bool ended = false;
	const char * name;
	bool debug;
public:

	TimeInterval(const char * aname, bool adebug) {
		name = aname;
		debug = adebug;
	}

	TimeInterval(std::string aname, bool adebug) {
		name = aname.c_str();
		debug = adebug;
	}

	TimeInterval(const char * aname) {
		name = aname;
		debug = debugtimers;
	}

	TimeInterval(std::string aname) {
		name = aname.c_str();
		debug = debugtimers;
	}


	void start() {
		if (started) {
			std::cout << "Timer " << name << " already started " << std::endl;
		}
		else {
			_start = clock();
			started = true;
		}
	}

	void end() {
		if (started) {
			_end = clock();
			ended = true;
		}
		else {
			std::cout << "Timer " << name << " not started " << std::endl;
		}
	}

	void print() {
		if (!ended) {
			std::cout << "Timer " << name << " not ended/started " << std::endl;
		}
		else {
			if (debug) {
				double full_elapsed_secs = double(_end - _start) / (CLOCKS_PER_SEC / 1000);
				std::cout << "Time(ms) " << name << ": " << full_elapsed_secs << std::endl;
			}
		}
	}



};

uint masks[32];
uint all = UINT32_MAX;
uint maxes[33];
uint none = 0;

#define ISBITON(value, ind) (value & masks[ind]) != 0

class element {

private:
	char count;
#define TURNON(ind) values = values | masks[ind]
#define TURNOFF(ind) values = values & ~masks[ind]
#define ISON(ind) (values & masks[ind]) != 0 
public:

	uint values;

	element(int size, bool initial) {
		if (initial) {
			values = all;
			count = size;
		}
		else {
			values = none;
			count = 0;
		}
	}

	element(int size, int onone) {
		values = none;
		count = 0;
		turnon(onone - 1);
		count = 1;
	}

	element(const element * el) {
		values = el->values;
		count = el->count;
	}

	bool iscertain() {
		return count == 1;
	}

	int getFirstElement() {
		for (int j = 0; j < 32; j++) {
			if (ISON(j)) {
				return j + 1;
			}
		}
		return -1;
	}

	int getcount() {
		return count;
	}

	bool turnoff(int index) {
		if (ISON(index)) {
			count--;
			TURNOFF(index);
		}
	}

	void turnon(int index) {
		if (!ISON(index)){
			count++;
			TURNON(index);
		}
	}

	void turnoff(element* other) {
		values = values & ~(other->values);
		count = __builtin_popcount(values);
	}

	void turnoffallbut(int one) {
		values = none;
		turnon(one);
		count = 1;
	}

	void turnon(element * other) {
		values = values | other->values;
		count = __builtin_popcount(values);
	}

	bool ison(int ind) {
		return ISON(ind);
	}

	string tostr(char zerochar, char size) {
		string buff;
		for (int i = 0; i < size; i++) {
			if (ISON(i)) {
				buff += (char)(i + (int)zerochar + 1);
				buff += ',';
			}
			else {
				buff += " ,";
			}
		}
		return buff;
	}

	
};

uint getOnlyOnes(element ** elements, int blocklength) {
	uint X = maxes[blocklength];
	uint Z = X;
	uint nA;
	for (int i = 0; i < blocklength; i++) {
		nA = ~elements[i]->values;
		Z = (Z & nA) | X;
		X =  (X & nA);
	}
	X = ~X&Z;
	return X;
}
//bit fo duplicaTION
uint getOnlyOnesNoCertains(vector<element**> rows, int dim, int blocklength) {
	uint X = maxes[blocklength];
	uint Z = X;
	uint nA;
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			if (!rows[i][j]->iscertain()) {//certians already removed from others
				nA = ~rows[i][j]->values;
				Z = (Z & nA) | X;
				X = (X & nA);
			}
		}
	}
	X = ~X&Z;
	return X;
}

class nrinRowColumnFilter {
public:
	void filter(element** row, int blocklength) {
		int count[blocklength];
		element* ref[blocklength];
		for (int i = 0; i < blocklength; i++) {
			count[i] = 0;
		}
		for (int i = 0; i < blocklength; i++) {
			for (int k = 0; k < blocklength; k++) {
				if (row[i]->ison(k)) {
					count[k]++;
					ref[k] = row[i];
				}
			}
		}
		for (int i = 0; i < blocklength; i++) {
			if (count[i] == 1 && !ref[i]->iscertain()) {
				ref[i]->turnoffallbut(i);
			}
		}

	}
};

class nrOfinABlock {
public:
	void filter(vector<element**> rows, int dim, int blocklength) {
		int count[blocklength];
		element* ref[blocklength];
		for (int i = 0; i < blocklength; i++) {
			count[i] = 0;
		}

		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				for (int k = 0; k < blocklength; k++) {
					if (rows[i][j]->ison(k)) {
						count[k]++;
						ref[k] = rows[i][j];
					}
				}
			}
		}
		for (int i = 0; i < blocklength; i++) {
			if (count[i] == 1 && !ref[i]->iscertain()) {
				ref[i]->turnoffallbut(i);
			}
		}
		



	}
};
class nrOfinABlockImporved {
public:
	void filter(vector<element**> rows, int dim, int blocklength) {
		uint onlyonce = getOnlyOnesNoCertains(rows, dim, blocklength);
	
		for (int k = 0; k < blocklength; k++) {
			if (ISBITON(onlyonce,k)) {
				for (int i = 0; i < dim; i++) {
					for (int j = 0; j < dim; j++) {
						if (rows[i][j]->ison(k)) {
							rows[i][j]->turnoffallbut(k);
						}
					}
				}
			}
		}

	}
};


class ruleOutOthers  {
public:
	void filter(vector<element**> rows, int dim, int blocklength) {
		element certain(blocklength,false);
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				if (rows[i][j]->iscertain()) {
					certain.turnon(rows[i][j]);
				}
			}
		}
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				if (!rows[i][j]->iscertain()) {
					 rows[i][j]->turnoff(&certain);
				}
			}
		}

	}
};

class BlockRowColumnFilter {
public:
	void filter(element** row, int blocklength) {
		element certain(blocklength, false);
		for (int i = 0; i < blocklength; i++) {
			if (row[i]->iscertain()) {
				certain.turnon(row[i]);
			}
		}
		for (int i = 0; i < blocklength; i++) {
			if (!row[i]->iscertain()) {
			 row[i]->turnoff(&certain);
			}
		}

	}
};

class excluderowblockcolumnfilter {
public:
	void filter(element*** row, int blocklength, int blockwidth) {
		vector<element*> markers;
		for (int i = 0; i < blocklength; i++) {
			markers.push_back(new element(blocklength, false));
		}
		for (int r = 0; r < blockwidth; r++) {
			for (int c = 0; c < blocklength; c++) {
				markers[r*blockwidth + (int)(c / blockwidth)]->turnon(row[r][c]);
			}
		}

		for (int c = 0; c < blockwidth; c++) {
			uint onlyonce;
			uint X = maxes[blocklength];
			uint Z = X;
			uint nA;
			for (int r = 0; r < blockwidth; r++) {
				nA = ~(markers[r*blockwidth + c]->values);
				Z = (Z & nA) | X;
				X = (X & nA);
			}
			onlyonce = ~X&Z;
			for (int i = 0; i < blocklength; i++) {
				if (ISBITON(onlyonce,i) == 1) {
					for (int r = 0; r < blockwidth; r++) {
						element * marker = markers[r*blockwidth + c];
						if (marker->ison(i)) {
							for (int c2 = 0; c2 < blocklength; c2++) {
								if ((int)(c2 / blockwidth) != c)//if not in the current block
									row[r][c2]->turnoff(i);
							}
							break;//it shoud happen only once
						}
					}
				}
			}
		}


	}
};

class Table {
private:
	vector<element*> cells;
	vector<vector<element**>> blocks;
	vector<element**> rows;
	vector<element**> columns;

	int length;
	int tablelength;
	int nrblocks;
	int blockheight;
	int blockwidth;
	bool waschange = false;
	

public:
	char zerochar;

	Table(int inputlength, char zerochararg) {
		length = inputlength;
		zerochar = zerochararg;
		tablelength = sqrt(length);
		blockheight = sqrt(tablelength);
		blockwidth = blockheight;
		nrblocks = length / tablelength; // :)
		blocks.resize(nrblocks);
		cells.resize(length);
		rows.resize(tablelength);
		columns.resize(tablelength);

	}

	Table(string input, char zerochararg) :Table((int)input.length(), zerochararg){
		
		for (int i = 0; i < length; i++) {
			if (input[i] != zerochar) {
				cells[i] = new element(tablelength, (int)input[i] - (int)zerochar);
			} else {
				cells[i] = new element(tablelength,true);
			}
		}
		initpointers();

	
	}

	Table(const Table *source) : Table(source->length, source->zerochar){
		for (int i = 0; i < length; i++) {
			cells[i] = new element(source->cells[i]);
		}
		//initpointers(); no need for teh if only doing recursion
	}

	void initpointers() {
		for (int i = 0; i < nrblocks; i++) {
			blocks[i].resize(blockheight);
			int blockrowstart = (i - (i % blockwidth)) * tablelength;
			for (int j = 0; j < blockheight; j++) {
				blocks[i][j] = cells.data() + ((i % blockwidth) * blockwidth + j * tablelength + blockrowstart);//pointer to each row of the block
			}
		}

		for (int i = 0; i < tablelength; i++) {
			rows[i] = cells.data() + i * tablelength;
		}

		for (int i = 0; i < tablelength; i++) {
			columns[i] = (element**)malloc(sizeof(element*) * tablelength);
			for (int j = 0; j < tablelength; j++) {
				columns[i][j] = *(cells.data() + j * tablelength + i);
			}
		}
	}

	void solve() {
		int loopnr = 0;
		int currcount = nrofremaining();
		int prevcount = 0;
		do {
			prevcount = currcount;
			applyBlockFilters();
			applyRowColumFilters();
			applyBlockRowColumnFilters();
			cout << loopnr << endl;
			loopnr++;
			currcount = nrofremaining();
		} while (currcount != prevcount);//if we ruled out something
	}

	void applyBlockFilters() {
		ruleOutOthers f;
#pragma omp parallel for 
		for (int j = 0; j < nrblocks; j++) {
			f.filter(blocks[j],blockwidth, tablelength);
		}
		nrOfinABlockImporved f2;
#pragma omp parallel for 
		for (int j = 0; j < nrblocks; j++) {
		    f2.filter(blocks[j], blockwidth, tablelength);
		}
	}

	void applyRowColumFilters() {
		BlockRowColumnFilter f;

#pragma omp parallel for 
		for (int i = 0; i < tablelength; i++) {
			 f.filter(rows[i], tablelength);
		}
#pragma omp parallel for 
		for (int i = 0; i < tablelength; i++) {
			 f.filter(columns[i], tablelength);
		}
		nrinRowColumnFilter f2;
#pragma omp parallel for 
		for (int i = 0; i < tablelength; i++) {
			 f2.filter(rows[i], tablelength);
		}
		bool waschangelocal4 = false;
#pragma omp parallel for 
		for (int i = 0; i < tablelength; i++) {
			 f2.filter(columns[i], tablelength);
		}

	}

	void applyBlockRowColumnFilters(){
		excluderowblockcolumnfilter f;
#pragma omp parallel for 
		for (int i = 0; i < blockwidth; i++) {
			f.filter(rows.data() + i * blockwidth, tablelength,blockwidth);
		}

#pragma omp parallel for 
		for (int i = 0; i < blockwidth; i++) {
			 f.filter(columns.data() + i * blockwidth, tablelength, blockwidth);
		}
	}

	void print() {
		for (int i = 0; i < tablelength; i++) {
			for (int j = 0; j < tablelength; j++) {
				cout << cells[i* tablelength + j]->tostr(zerochar, tablelength) << "\t";
			}
			cout << endl;
		}
		cout << "nrpos " << nrofpossibilities()  <<  " nrofsures " << nrofsures() << endl;
	}

	long double nrofpossibilities() {
		long double nr = 1;
		for (int i = 0; i < tablelength; i++) {
			for (int j = 0; j < tablelength; j++) {
				nr *= cells[i* tablelength + j]->getcount();
			}
		}
		return nr;
	}

	string tosureonesstring() {
		string buff;
		for (int i = 0; i < length; i++) {
			if (cells[i]->iscertain()) {
				buff += (char)(cells[i]->getFirstElement() + (int)zerochar);
			}
			else {
				buff += zerochar;
			}
		}
		return buff;
	}

	int nrofsures() {
		int sum = 0;
		for (int i = 0; i < length; i++) {
			if (cells[i]->iscertain()) {
				sum++;
			}
		}
		return sum;
	}

	int nrofremaining() {
		int sum = 0;
		for (int i = 0; i < length; i++) {
			sum += cells[i]->getcount();
		}		
		return sum;
	}
/*	bool isSolved() {
		for (int i = 0; i < length; i++) {
			if (!cells[i]->iscertain()) {
				return false;
			}
		}
		return true;
	}

	bool solveBackTrack()
	{
		// Készen vagyunk?
		if (isSolved())
		{
			return true;
		}

		// Keressünk egy pozíciót, amely még nincs kitöltve
		for (int y = 0; y < tablelength; ++y)
		{
			for (int x = 0; x < tablelength; ++x)
			{
				// Nincs még kitöltve?
				if (!cells[y * tablelength + x]->iscertain())
				{
					// Keressünk egy értéket, amely megfelel a szabályoknak
					for (int n = 0; n < blockwidth; ++n)
					{
						// Beírható az adott pozícióba?
						if (isAllowed(n, x, y))
						{
							// Másoljuk le a táblát
							Table tmpTable(this);
							// Írjuk bele az új értéket
							tmpTable.set(n, x, y);
							// Próbáljuk megoldani az új táblát
							if (tmpTable.solveBackTrack())
							{
								// Megoldás
								*this = tmpTable;
								return true;
							}
						}
					}
				}
				// Nem tudtunk értéket írni a cellába, így lépjünk vissza
				if (!cells[y * tablelength + x]->iscertain()) return false;
			}
		}
	
		return false;
	}

	void set(char val, int x, int y)
	{
		cells[y * tablelength + x]->turnoffallbut(val);
	}


	bool isAllowed(char val, int x, int y)
	{
		if (!cells[y * tablelength + x]->ison(val)) {
			return false;
		}
		// Azonos sorban vagy oszlopban csak egy 'val' lehet
		for (int i = 0; i < tablelength; ++i)
		{
			if (cells[y * tablelength + i]->ison(val) && cells[y * tablelength + i]->iscertain()) return false;
			if (cells[i * tablelength + x]->ison(val) && cells[i * tablelength + x]->iscertain()) return false;
		}

		// Az adott blockxblock-as cellában csak egy 'val' lehet
		int cellBaseX = blockwidth * (int)(x / blockwidth);
		int cellBaseY = blockwidth * (int)(y / blockwidth);
		for (int y = cellBaseY; y < cellBaseY + blockwidth; ++y)
		{
			for (int x = cellBaseX; x < cellBaseX + blockwidth; ++x)
			{
				if (cells[y * tablelength + x]->ison(val) && cells[y * tablelength + x]->iscertain()) return false;
			}
		}

		return true;
	}*/
};


void testof() {
	uint tests[8][5] = { 
		{1,1,1,1,0},//0
		{1,0,0,1,0},//1
		{0,0,0,0,0},//2
		{0,0,1,0,1},//3
		{1,0,0,0,1},//4
		{0,0,0,1,1},//5
		{0,1,0,0,1},//6
		{0,0,1,1,0},//7 
	};

	element ** els = new element*[4];
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 4; j++) {
			els[j] = new element(4, false);
			if (tests[i][j] == 1) {
				els[j]->turnon(0);
			}
		}
		if(getOnlyOnes(els, 4) != tests[i][4]){
			cout << i << " not correct" << endl;
		}
	}

}


int main()
{

	for (int i = 0; i < 32; i++) {
		masks[i] = 1 << i;
	}
	maxes[0] = 0;
	maxes[1] = 1;
	for (int i = 2; i <= 32; i++) {
		maxes[i] = maxes[i-1] * 2;
	}
	for (int i = 2; i <= 32; i++) {
		maxes[i] --;
	}
	/*testof();
	return 0;*/
	
	//Table t("@BM@EO@@@P@@C@@@@@IAB@@@O@GFK@@@@GCP@AFN@EKIB@J@FEKNIL@GBJC@OP@@@@FCJN@@@HLOGE@@@@@@@DKPG@@@@O@@@OP@@@@E@@BKA@MD@@@GO@I@@@@@@CK@@NG@HIL@POA@@@@@ALJ@@@@BIC@@@KO@@@OI@E@@@@DB@@P@@PB@@@@@JL@@I@@@@CNBAM@I@@OLP@@@@J@KGB@LNA@@@IFOE@@@@JC@@@IP@GBL@IL@@KE@@@@GN@AC", (char)64);
	//Table t("000000012500008000000700000600120000700000450000030000030000800000500700020000000",'0');
	Table t("000801000000000043500000000000070800020030000000000100600000075003400000000200600",'0');
	//Table t("000000600079624008564800001605090783000003000823040010230060007001200060790050002",'0');
	//t.print();
	TimeInterval solvetfull("solvefull");
	TimeInterval solvet("solve");
	solvet.start();
	solvetfull.start();
	t.solve();
	solvet.end();
	TimeInterval solvet2("solvebr");
	solvet2.start();
	cout << t.tosureonesstring().c_str() << endl;
	Solver solver2(t.tosureonesstring(),t.zerochar);
	t.solveBackTrack();
	solvetfull.end();
	solvet2.end();
	solvet2.print();
	t.print();
	solvet.print();
	solvetfull.print();
	t.print();
	//solver2.print(std::cout);


	//Solver solver("000000012500008000000700000600120000700000450000030000030000800000500700020000000");
	Solver solver("@BM@EO@@@P@@C@@@@@IAB@@@O@G@K@@@@GCP@AFN@EKIB@J@FEKNIL@GBJC@OP@@@@FCJN@@@HL@@E@@@@@@@D@@G@@@@O@@@OP@@@@E@@BKA@MD@@@GO@I@@@@@@CK@@NG@H@L@POA@@@@@ALJ@@@@B@C@@@KO@@@OI@E@@@@DB@@P@@PB@@@@@JL@@I@@@@CN@@M@I@@OLP@@@@J@K@B@LNA@@@IFOE@@@@JC@@@IP@GBL@IL@@KE@@@@GN@AC", (char)64);

	//Solver solver("000801000000000043500000000000070800020030000000000100600000075003400000000200600",'0');
	std::cout << "Problem:" << std::endl << std::endl;
	solver.print(std::cout);
	std::cout << std::endl << "-----------------------------------------" << std::endl;
	std::cout << "Solution:" << std::endl << std::endl;;
	TimeInterval solvebackt("solvebacktarck");
	solvebackt.start();
//	solver.solveBackTrack();
	solvebackt.end();
	solver.print(std::cout);
	solvebackt.print();

	cout.flush();
	usleep(2000);

    return 0;
}