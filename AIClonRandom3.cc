#include "Player.hh"
#include <vector>
#include <map>
#include <stack>

using namespace std;

#define PLAYER_NAME ClonRandom3

struct PLAYER_NAME : public Player {

	#define MIN_MISSILES 4
	#define STRICT_SIMULATION true
	#define BIG_SIMULATION true
	#define SIMULATION_ON_SCAN true
	#define LIMIT_POS 0.9 //[0-1] 1 = no limit
	#define MAX_LEVEL_BFS 20
	#define N_TARGETS_BFS 2

	static Player* factory () {
		return new PLAYER_NAME;
	}

	int m_universe;

	const Dir INVALID_DIR = {-1, -1};

	struct Data_pos {
		Pos prev;
		int missiles;
	};

	typedef pair<Pos, Data_pos> Nodo;

	struct compare {
		bool operator() (const Pos& p1, const Pos& p2) const{
			if(first(p1) != first(p2)) {
				return first(p1) < first(p2);
			} else {
				return second(p1) < second(p2);
			}
		}
	};

	typedef map<Pos, Data_pos, compare> Nodos;

	enum TARGET_TYPE{
		POINTS,
		MISSILES,
		ANY,
		ENEMY,
		REPOSITION
	};

	struct Target {
		Nodo n;
		stack<Pos> route;
		int missiles_nec;
		int rounds_nec;
		bool ok;
		TARGET_TYPE type;
		Pos next_pos;

		void set_route(Nodos &nodos, const Nodo &n) {
			if(n.first != n.second.prev) {
				route.push(n.first);
				Nodos::iterator i = nodos.find(n.second.prev);
				if(i != nodos.end()) {
					Nodo aux = *i;
					nodos.erase(i);
					set_route(nodos, aux);
				}
			}
		}

		void add_step(const Pos &p) {
			route.push(p);
		}

		void reset_route() {
			 route = stack<Pos>();
		}
	};

	class Targets {
	public:
		Targets() {}

		Targets(int size, int offset, int m_universe) {
			this->m_universe = m_universe;
			this->offset = offset;
			Target t;
			t.ok = false;
			v = vector<Target>(size, t);
		}

		Target& operator[](Starship_Id id) {
			return v[id-offset];
		}

		bool avaliable(Starship_Id id, const Pos &p) {
			for(int i = 0; i < (int)v.size(); ++i) {
				if(id-offset != i && v[i].ok && first(v[i].n.first) == first(p)%m_universe && second(v[i].n.first)%m_universe == second(p)%m_universe) {
					return false;
				}
			}
			return true;
		}
		int offset = 0;
	private:
		int m_universe;
		vector<Target> v;
	};

	Targets targets;

	typedef vector<vector<Cell> > Stage;

	//TODO add const methods
	class Simulation {
	public:
		Stage stage;

		Simulation() {}

		Simulation(Pos from, int n_universe, int m_universe, bool strict, bool big, int border) {
			this->m_universe = m_universe;
			this->strict = strict;
			reduce_pos(from);
			int n = 3+2*big;
			int m = 5+border;
			int i = first(from)-1-big;
			int j = second(from)-4;
			if(i < 0) {
				n = n+i;
				i = 0;
			}
			ref = {i, j};
			reduce_pos(ref);
			while(i+n > n_universe) {
				--n;
			}
			stage = Stage(n, vector<Cell>(m));
		}

		void run(const vector<Dir> &all_dirs) {
			for(int i = 0; i < (int)stage.size(); ++i) {
				for(int j = (int)stage[0].size()-1; j >= 0; --j) {
					if(stage[i][j].type == MISSILE) {
						stage[i][j].type = EMPTY;
						bool crash = false;
						for(int k = 1; !crash && k <= 2; ++k) {
							if(j+k < (int)stage[0].size() && stage[i][j+k].type != EMPTY && stage[i][j+k].sid < 0) {
								crash = true;
								stage[i][j+k].type = EMPTY;
							}
						}
						if(!crash && j+2 < (int)stage[0].size()) {
							stage[i][j+2].type = MISSILE;
						}
					}
				}
			}
			if(strict) {
				for(int j = (int)stage[0].size()-1; j >= 0; --j) {
					for(int i = 0; i < (int)stage.size(); ++i) {
						if(stage[i][j].type == STARSHIP && stage[i][j].sid == -1) {
							stack<pair<pair<int, int>, Cell> > cells = stack<pair<pair<int, int>, Cell> >();
							for(int d = 0; d < (int)all_dirs.size(); ++d) {
								pair<int, int> ij = {i, j};
								if(check_movement(ij, all_dirs[d], stage[i][j].mid > 0)) {
									ij.first += first(all_dirs[d]);
									ij.second += second(all_dirs[d]);
									if(affects(ij)) {
										pair<pair<int, int>, Cell> new_pos;
										new_pos.first = ij;
										if(new_pos.second.type == ASTEROID) {
											--new_pos.second.mid;
										}
										new_pos.second = stage[i][j];
										new_pos.second.sid = -2,
										cells.push(new_pos);
									}
								}
							}
							while(!cells.empty()) {
								stage[cells.top().first.first][cells.top().first.second] = cells.top().second;
								cells.pop();
							}
						}
					}
				}
			}
		}

		bool can_move(Pos p, Dir d, bool has_missiles) const {
			reduce_pos(p);
			pair<int, int> ij = get_ij(p);
			if(!check_movement(ij, d, has_missiles && stage[get_ij(p+d).first][get_ij(p+d).second].sid != -2)) {
				return false;
			}
			p = p+d;
			for(int j = -2; j < 0; ++j) {
				d = {0, j};
				pair<int, int> ij = get_ij(p+d);
				if(affects(ij)) {
					Cell c = stage[ij.first][ij.second];
					if(c.type == MISSILE || (c.type == STARSHIP && (c.sid == -1 || (c.sid == -2 && strict)) && c.mid > 0)) {
						if(j == -1 || stage[ij.first][ij.second+1].type == EMPTY || stage[ij.first][ij.second+1].type == STARSHIP) {
							return false;
						}
					}
				}
			}
			return true;
		}

		Dir get_free_dir(Pos from, const vector<Dir> &all_dirs, bool slow, bool has_missiles) {
			reduce_pos(from);
			for(int i = 0; i < (int)all_dirs.size(); ++i) {
				int d = priorize_dir(all_dirs, i, slow);
				pair<int, int> ij = get_ij(from+all_dirs[d]);
				if(affects(ij) && can_move(from, all_dirs[d], has_missiles && stage[ij.first][ij.second].sid != -2)) {
					return all_dirs[d];
				}
			}
			if(strict) {
				//cerr << "no hay dir libre, se baja el nivel de strict" << endl;
				strict = false;
				Dir d = get_free_dir(from, all_dirs, slow, has_missiles);
				strict = true;
				return d;
			} else {
				//cerr << "no hay ninguna dir libre, default" << endl;
				return DEFAULT;
			}
		}

		void add_starship(Starship_Id id, Pos p) {
			reduce_pos(p);
			pair<int, int> ij = get_ij(p);
			if(affects(ij)) {
				stage[ij.first][ij.second].type = STARSHIP;
				stage[ij.first][ij.second].sid = id;
			}
		}

		void set_cell(Pos p, const Cell &c) {
			reduce_pos(p);
			pair<int, int> ij = get_ij(p);
			if(affects(ij)) {
				stage[ij.first][ij.second] = c;
			}
		}

		Cell get_cell(Pos p) {
			reduce_pos(p);
			pair<int, int> ij = get_ij(p);
			if(affects(ij)) {
				return stage[ij.first][ij.second];
			}
			Cell c;
			c.type = EMPTY;
			return c;
		}

		Pos ref;

	private:
		int m_universe;
		bool strict;

		void reduce_pos(Pos &p) const {
			int j = second(p)%m_universe;
			if(j < 0) {
				j = m_universe+j;
			}
			p = {first(p), j};
		}

		bool affects(const pair<int, int> &ij) const {
			return ij.first >= 0 && ij.first < (int)stage.size() && ij.second >= 0 && ij.second < (int)stage[0].size();
		}

		pair<int, int> get_ij(Pos p) const {
			p = p-ref;
			reduce_pos(p);
			pair<int, int> ij;
			ij.first = first(p);
			ij.second = second(p);
			return ij;
		}

		bool cell_ok(pair<int, int> ij, const Dir &d) const {
			ij.first += first(d);
			ij.second += second(d);
			if(affects(ij)) {
				Cell c = stage[ij.first][ij.second];
				return c.type != ASTEROID && c.type != MISSILE && (c.type != STARSHIP || (!strict && c.sid == -2));
			} else {
				return false;
			}
		}

		bool check_movement(const pair<int, int> &ij, const Dir &d, bool has_missiles) const {
			if(d == DEFAULT) {
				return cell_ok(ij, d) || (has_missiles && stage[ij.first][ij.second+1].sid < 0);
			} else if(!cell_ok(ij, d)) {
				return false;
			} else if(d == FAST) {
				return cell_ok(ij, DEFAULT);
			} else if(d == FAST_UP) {
				return check_movement(ij, UP, has_missiles) && check_movement(ij, FAST, has_missiles);
			} else if(d == FAST_DOWN) {
				return check_movement(ij, DOWN, has_missiles) && check_movement(ij, FAST, has_missiles);
			} else if(d == UP) {
				return cell_ok(ij, SLOW_UP) && cell_ok(ij, DEFAULT);
			} else if(d == DOWN) {
				return cell_ok(ij, SLOW_DOWN) && cell_ok(ij, DEFAULT);
			} else {
				return true;
			}
		}
	};

	Simulation simulation;

	void load_simulation(const Starship &s) {
		int border = 2;
		Dir aux = {0, 2};
		if(!within_window(s.pos+aux, round())) {
			--border;
			aux = {0, 1};
			if(!within_window(s.pos+aux, round())) {
				--border;
			}
		}
		simulation = Simulation(s.pos, number_rows(), m_universe, STRICT_SIMULATION, BIG_SIMULATION, border);
		Dir d;
		for(int i = 0; i < (int)simulation.stage.size(); ++i) {
			for(int j = 0; j < (int)simulation.stage[0].size(); ++j) {
				d = {i, j};
				simulation.stage[i][j] = cell(simulation.ref+d);
				if(simulation.stage[i][j].type == STARSHIP) {
					if(player_of(simulation.stage[i][j].sid) != me()) {
						simulation.stage[i][j].mid = starship(simulation.stage[i][j].sid).nb_miss;
						simulation.stage[i][j].sid = -1;
					}
				}
			}
		}
		for(Starship_Id id_aux = begin(me()); id_aux != end(me()); ++id_aux) {
			if(id_aux != s.sid && targets[id_aux].ok) {
				simulation.add_starship(id_aux, targets[id_aux].next_pos);
			}
		}
		//print_simulation(simulation);
		simulation.run(all_dirs);
		//print_simulation(simulation);
	}

	vector<Dir> all_dirs;

	bool cell_ok(const Pos &p) {
		return cell(p).type != ASTEROID;
	}

	bool check_movement(const Pos &p, const Dir &d, bool has_missiles) {
		if(d == DEFAULT) {
			return cell_ok(p+d) || has_missiles;
		} else if(!cell_ok(p+d)) {
			return false;
		} else if(d == FAST) {
			return cell_ok(p+DEFAULT);
		} else if(d == FAST_UP) {
			return check_movement(p, UP, has_missiles) && check_movement(p, FAST, has_missiles);
		} else if(d == FAST_DOWN) {
			return check_movement(p, DOWN, has_missiles) && check_movement(p, FAST, has_missiles);
		} else if(d == UP) {
			return cell_ok(p+SLOW_UP) && cell_ok(p+DEFAULT);
		} else if(d == DOWN) {
			return cell_ok(p+SLOW_DOWN) && cell_ok(p+DEFAULT);
		} else {
			return true;
		}
	}

	Target choose_target(const vector<Target> &v, int size) {
		//TODO escoger el target que tenga otro target del mismo tipo en menos de n (peque�a) num de rondas desde ahi.
		int k = -1;
		int min_j = 0;
		for(int i = 0; i < size; ++i) {
			if(k == -1 || second(v[i].n.first) < min_j) {
				min_j = second(v[i].n.first);
				k = i;
			}
		}
		return v[k];
	}

	bool is_target(const Target &t, CType type) {
		if(t.type == REPOSITION && (type == POINT_BONUS || type == MISSILE_BONUS)) {
			return true;
		} else if(t.type == POINTS && type == POINT_BONUS) {
			return true;
		} else if(t.type == MISSILES && type == MISSILE_BONUS) {
			return true;
		} else if(t.type == ANY && (type == MISSILE_BONUS || type == POINT_BONUS)) {
			return true;
		} else if(t.type == ENEMY && type == STARSHIP) {
			return true;
		}
		return false;
	}

	static int priorize_dir(const vector<Dir> &all_dirs, int d, bool slow) {
		if(!slow) {
			return d;
		} else {
			return all_dirs.size()-1-d;
		}
	}

	//TODO realizar adaptacion alternativa algoritmo dijkstra: busqueda+camino minimo a la vez
	bool scan_target(const Starship &s, Target &t) {
		t.reset_route();
		bool slow = t.type == REPOSITION;
		vector<Target> v = vector<Target>(N_TARGETS_BFS);
		Nodos visited = Nodos();
		queue<Nodo> q;
		Nodo n;
		n.first = s.pos;
		n.second.missiles = s.nb_miss;
		n.second.prev = n.first;
		q.push(n);
		visited.insert(n);
		int r = 0;
		int i = 0;
		int current_level = 1;
		int next_level = 0;
		Nodo act;
		while(!q.empty() && i < N_TARGETS_BFS && r <= MAX_LEVEL_BFS) {
			act = q.front();
			if(is_target(t, cell(act.first).type)) {
				if(targets.avaliable(s.sid, act.first)) {
					v[i].n = act;
					v[i].rounds_nec = r;
					v[i].missiles_nec = s.nb_miss-act.second.missiles;
					++i;
				}
			}
			for(int j = 0; j < (int)all_dirs.size(); ++j) {
				int d = priorize_dir(all_dirs, j, slow);
				Pos pos_new = act.first+all_dirs[d];
				if((visited.find(pos_new) == visited.end() || pos_new == act.first)
						&& within_window(pos_new, round()+r)
						&& check_movement(act.first, all_dirs[d], act.second.missiles > 0)) {
					if(!SIMULATION_ON_SCAN || r != 0 || simulation.can_move(act.first, all_dirs[d], s.nb_miss > 0)) {
						n.first = pos_new;
						n.second.prev = act.first;
						n.second.missiles = act.second.missiles - !cell_ok(pos_new);
						visited.insert(n);
						q.push(n);
						++next_level;
					}
				}
			}
			q.pop();
			--current_level;
			if(current_level == 0) {
				current_level = next_level;
				next_level = 0;
				++r;
			}
		}
		if(i > 0) {
			t = choose_target(v, i);
			t.set_route(visited, t.n);
			return true;
		} else {
			return false;
		}
	}

	void print_pos(const Pos &p) {
		//cerr << '(' << first(p) << ',' << second(p)%m_universe << ')' << endl;
	}

	void print_route(const Starship &s) {
		stack<Pos> pila = targets[s.sid].route;
		while(!pila.empty()) {
			//print_pos(pila.top());
			pila.pop();
		}
	}

	void print_simulation(const Simulation &s) {
		for(int i = 0; i < (int)s.stage.size(); ++i) {
			for(int j = 0; j < (int)s.stage[0].size(); ++j) {
				Cell c = s.stage[i][j];
				if(c.type == ASTEROID) {
					//cerr << 'X';
				} else if(c.type == POINT_BONUS) {
					//cerr << 'P';
				} else if(c.type == MISSILE_BONUS) {
					//cerr << 'B';
				} else if(c.type == MISSILE) {
					//cerr << 'M';
				} else if(c.type == STARSHIP) {
					if(c.sid < 0) {
						//cerr << 'S';
					} else {
						//cerr << c.sid;
					}
				} else {
					//cerr << '.';
				}
			}
			//cerr << endl;
		}
		//cerr << endl;
	}

	//TODO estructura de datos para almacenar naves enemigas
	//TODO posible campo rondas espera sin adjudicar target para evitar hacer un bfs cada ronda si en la anterior no habia nada

	Dir get_dir(Pos from, Pos to) {
		from = {first(from), second(from)%m_universe};
		to = {first(to), second(to)%m_universe};
		Dir d = to-from;
		if(second(d) < 0) {
			d = {first(d), m_universe+second(d)};
		}
		for(int i = 0; i < (int)all_dirs.size(); ++i) {
			if(d == all_dirs[i]) {
				return d;
			}
		}
		return INVALID_DIR;
	}

	bool dispara(const Starship &s) {
		Dir d;
		for(int j = 1; j <= 2; ++j) {
			d = {0, j};
			if(cell(s.pos+d).type == STARSHIP && player_of(cell(s.pos+d).sid) == me()) {
				if(j == 1) {
					targets[s.sid].next_pos = s.pos+SLOW;
				} else {
					targets[s.sid].next_pos = s.pos+DEFAULT;
				}
				return false;
			} else if(cell(s.pos+d).type != EMPTY) {
				shoot(s.sid);
				targets[s.sid].next_pos = s.pos+DEFAULT;
				return true;
			}
		}
		shoot(s.sid);
		targets[s.sid].next_pos = s.pos+DEFAULT;
		return true;
	}

	void recalculate_route(const Starship &s) {
		targets[s.sid].type = ANY;
		targets[s.sid].ok = scan_target(s, targets[s.sid]);
		if(!targets[s.sid].ok) {
			targets[s.sid].reset_route();
			Dir d = simulation.get_free_dir(s.pos, all_dirs, targets[s.sid].type == REPOSITION, s.nb_miss > 0);
			targets[s.sid].route.push(s.pos+d);
		}
	}

	void check_safe(const Starship &s) {
		Dir d = get_dir(s.pos, targets[s.sid].route.top());
		if(!simulation.can_move(s.pos, d, s.nb_miss > 0)) {
			//cerr << "Movimiento no seguro, recalcular ruta" << endl;
			recalculate_route(s);
		} else {
			//cerr << "Movimiento previsto seguro, no hay cambios" << endl;
		}
	}

	void refresh_target(const Starship &s) {
		if(!targets[s.sid].ok || !is_target(targets[s.sid], cell(targets[s.sid].n.first).type) || get_dir(s.pos, targets[s.sid].route.top()) == INVALID_DIR) {
			TARGET_TYPE type = POINTS;
			if(s.nb_miss < MIN_MISSILES) {
				type = MISSILES;
				//cerr << "tiene que buscar misiles" << endl;
			}
			if(second(s.pos)%m_universe >= round()%m_universe+number_window_columns()*LIMIT_POS) {
				type = REPOSITION;
				//cerr << "Demasiado adelante" << endl;
			}
			targets[s.sid].type = type;
			targets[s.sid].ok = scan_target(s, targets[s.sid]);
			if(targets[s.sid].ok == true) {
				//cerr << "Se le adjudica objetivo: " << targets[s.sid].type << " en ";
				//print_pos(targets[s.sid].n.first);
			}
		}
		if(!targets[s.sid].ok) {
			targets[s.sid].route.push(s.pos+DEFAULT);
			targets[s.sid].ok = true;
			//cerr << "No se le ha encontrdo objetivo, se le a�ade dir default" << endl;
		}
	}

	///
	/// \brief play    asigna un objetivo a cada nave al principio y mueve a cada nave seg?n su objetivo. Tambi�n actualiza
	/// la posici�n de las naves enemigas en funci�n de sus antiguas posiciones para no tener que hacer una b?squeda en cada turno.
	///

	virtual void play() {
		//cerr << "---------------------------" << endl;
		//cerr << "RONDA " << round() << endl;
		if(round() == 0) {
			m_universe = number_universe_columns();
			all_dirs = {FAST, FAST_UP, FAST_DOWN, DEFAULT, UP, DOWN, SLOW_UP, SLOW_DOWN, SLOW};
			targets = Targets(number_starships_per_player(), begin(me()), m_universe);
		}
		for(Starship_Id id = begin(me()); id != end(me()); ++id) {
			const Starship s = starship(id);
			if(s.alive) {
				//cerr << "Turno de starship " << id << endl;
				//cerr << "Pos: ";
				//print_pos(s.pos);
				if(targets[id].route.empty()) {
					//cerr << "No tiene ruta" << endl;
					targets[id].ok = false;
				}

				load_simulation(s);
				refresh_target(s);
				check_safe(s);

				//cerr << "Su ruta actual es: " << endl;
				//print_route(s);

				if(false) {
				} else if (cell(targets[id].route.top()).type == ASTEROID || cell(targets[id].route.top()).type == STARSHIP) {
					dispara(s);
				} else {
					move(id, get_dir(s.pos, targets[id].route.top()));
				}
				targets[id].next_pos = targets[id].route.top();
				targets[id].route.pop();
			} else {
				//while(true);
				targets[id].ok = false;
			}
		}
	}
};

RegisterPlayer(PLAYER_NAME);

