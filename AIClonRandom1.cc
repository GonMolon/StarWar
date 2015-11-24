#include "Player.hh"
#include <vector>
#include <map>
#include <stack>

using namespace std;

#define PLAYER_NAME ClonRandom1

struct PLAYER_NAME : public Player {

	static Player* factory () {
		return new PLAYER_NAME;
	}

	int m_universe;

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
		POSITION
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
	};

	class Targets {
	public:
		Targets() {}

		Targets(int size, int offset, int m) {
			this->m = m;
			this->offset = offset;
			Target t;
			t.ok = false;
			v = vector<Target>(size, t);
		}

		Target& operator[](Starship_Id id) {
			return v[id-offset];
		}

		bool avaliable(Starship_Id id, const Pos &p) {
			for(int i = 0; i < v.size(); ++i) {
				if(id-offset != i && v[i].ok && first(v[i].n.first) == first(p) && second(v[i].n.first)%m == second(p)) {
					return false;
				}
			}
			return true;
		}
		int offset = 0;
	private:
		int m;
		vector<Target> v;
	};

	Targets targets;

	typedef vector<vector<Cell> > Stage;

	//TODO add const methods
	class Simulation {
	public:
		Stage stage;

		Simulation() {}

		Simulation(Pos from, Pos to, int n_universe, int m_universe, int border) {
			this->m_universe = m_universe;
			reduce_pos(from);
			reduce_pos(to);
			int i = first(from)-1;
			int j = second(to)-4;
			if(i < 0) {
				i = 0;
			}
			ref = {i, j};
			reduce_pos(ref);
			int n, m;
			if(first(from) == n_universe-1 || first(from) == 0) {
				n = 2;
			} else {
				n = 3;
			}
			Pos aux = from-ref;
			reduce_pos(aux);
			m = second(aux)+border+1;
			stage = Stage(n, vector<Cell>(m));
		}

		void run(bool strict, const vector<Dir> &all_dirs) {
			for(int i = 0; i < stage.size(); ++i) {
				for(int j = stage[0].size()-1; j >= 0; --j) {
					if(stage[i][j].type == MISSILE) {
						stage[i][j].type = EMPTY;
						bool crash = false;
						for(int k = 1; !crash && k <= 2; ++k) {
							if(j+k < stage[0].size() && stage[i][j+k].type != EMPTY) {
								crash = true;
								stage[i][j+k].type = EMPTY;
							}
						}
						if(!crash && j+2 < stage[0].size()) {
							stage[i][j+2].type = MISSILE;
						}
					}
				}
			}
			if(strict) {
				for(int j = stage[0].size()-1; j >= 0; --j) {
					for(int i = 0; i < stage.size(); ++i) {
						if(stage[i][j].type == STARSHIP && stage[i][j].sid == -1) {
							for(int d = 0; d < all_dirs.size(); ++d) {
								pair<int, int> ij;
								ij.first = i+first(all_dirs[d]);
								ij.second = j+second(all_dirs[d]);
								if(affects(ij) && stage[ij.first][ij.second].type == EMPTY) {
									stage[ij.first][ij.second] = stage[i][j];
								}
							}
						}
					}
				}
			}
		}

		bool can_move(Pos p, Dir d) const {
			reduce_pos(p);
			if(!check_movement(p, d)) {
				return false;
			}
			p = p+d;
			for(int j = -2; j < 0; ++j) {
				d = {0, j};
				pair<int, int> ij = get_ij(p+d);
				if(affects(ij)) {
					Cell c = stage[ij.first][ij.second];
					if(c.type == MISSILE || (c.type == STARSHIP && c.sid == -1)) {
						return false;
					}
				}
			}
			return true;
		}

		Dir get_free_dir(Pos from, const vector<Dir> &all_dirs) const {
			reduce_pos(from);
			for(int d = 0; d < all_dirs.size(); ++d) {
				pair<int, int> ij = get_ij(from+all_dirs[d]);
				if(affects(ij) && can_move(from, all_dirs[d])) {
					return all_dirs[d];
				}
			}
			return DEFAULT;
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

		void reduce_pos(Pos &p) const {
			int j = second(p)%m_universe;
			if(j < 0) {
				j = m_universe+j;
			}
			p = {first(p), j};
		}

		bool affects(const pair<int, int> &ij) const {
			return ij.first >= 0 && ij.first < stage.size() && ij.second >= 0 && ij.second < stage[0].size();
		}

		pair<int, int> get_ij(Pos p) const {
			p = p-ref;
			reduce_pos(p);
			pair<int, int> ij;
			ij.first = first(p);
			ij.second = second(p);
			return ij;
		}

		bool cell_ok(const Pos &p) const {
			pair<int, int> ij = get_ij(p);
			if(affects(ij)) {
				Cell c = stage[ij.first][ij.second];
				return c.type != ASTEROID && c.type != MISSILE && c.type != STARSHIP;
			} else {
				return false;
			}
		}

		bool check_movement(const Pos &p, const Dir &d) const {
			if(!cell_ok(p+d)) {
				return false;
			} else if(d == FAST) {
				return cell_ok(p+DEFAULT);
			} else if(d == FAST_UP) {
				return check_movement(p, UP) && check_movement(p, FAST);
			} else if(d == FAST_DOWN) {
				return check_movement(p, DOWN) && check_movement(p, FAST);
			} else if(d == UP) {
				return cell_ok(p+SLOW_UP) && cell_ok(p+DEFAULT);
			} else if(d == DOWN) {
				return cell_ok(p+SLOW_DOWN) && cell_ok(p+DEFAULT);
			} else {
				return true;
			}
		}
	};

	void load_simulation(Starship_Id id, Simulation &simulation) {
		Dir d;
		for(int i = 0; i < simulation.stage.size(); ++i) {
			for(int j = 0; j < simulation.stage[0].size(); ++j) {
				d = {i, j};
				simulation.stage[i][j] = cell(simulation.ref+d);
				if(simulation.stage[i][j].type == STARSHIP) {
					if(player_of(simulation.stage[i][j].sid) != me()) {
						simulation.stage[i][j].sid = -1;
					}
				}
			}
		}
		for(Starship_Id id_aux = begin(me()); id_aux != end(me()); ++id_aux) {
			if(id_aux != id && targets[id_aux].ok) {
				simulation.add_starship(id_aux, targets[id_aux].next_pos);
			}
		}
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

	///
	/// \brief choose_target sirve para escoger un objetivo y un recorrido para
	/// \param s
	/// \return devuelve el objetivo seleccionado
	///

	Target choose_target(const vector<Target> &v, int size) {
		//TODO implementar choose_target
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

	bool is_target(const Target &t, const Pos &p) {
		if(t.type == POSITION && t.n.first == p) {
			return true;
		} else if(t.type == POINTS && cell(p).type == POINT_BONUS) {
			return true;
		} else if(t.type == MISSILES && cell(p).type == MISSILE_BONUS) {
			return true;
		} else if(t.type == ANY && (cell(p).type == MISSILE_BONUS || cell(p).type == POINT_BONUS)) {
			return true;
		} else if(t.type == ENEMY && cell(p).type == STARSHIP) {
			return true;
		}
		return false;
	}

	#define MAX_LEVEL 15
	#define N_TARGETS 2

	//TODO realizar adaptacion alternativa algoritmo dijkstra: busqueda+camino minimo a la vez
	bool scan_target(const Starship &s, Target &t) {
		t.route = stack<Pos>();
		vector<Target> v = vector<Target>(N_TARGETS);
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
		while(!q.empty() && i < N_TARGETS) { //&& r <= MAX_LEVEL) {
			act = q.front();
			if(is_target(t, act.first)) {
				if(targets.avaliable(s.sid, act.first)) {
					v[i].n = act;
					v[i].rounds_nec = r;
					v[i].missiles_nec = s.nb_miss-act.second.missiles;
					++i;
				}
			}
			for(int j = 0; j < all_dirs.size(); ++j) {
				Pos pos_new = act.first+all_dirs[j];
				if(visited.find(pos_new) == visited.end()
						&& within_window(pos_new, round()+r)
						&& check_movement(act.first, all_dirs[j], act.second.missiles > 0)) {
					n.first = pos_new;
					n.second.prev = act.first;
					n.second.missiles = act.second.missiles - !cell_ok(pos_new);
					visited.insert(n);
					q.push(n);
					++next_level;
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
		//cerr << '(' << first(p) << ',' << second(p)%m << ')' << endl;
	}

	void print_route(const Starship &s) {
		stack<Pos> pila = targets[s.sid].route;
		while(!pila.empty()) {
			//print_pos(pila.top());
			pila.pop();
		}
	}

	void print_simulation(const Simulation &s) {
		for(int i = 0; i < s.stage.size(); ++i) {
			for(int j = 0; j < s.stage[0].size(); ++j) {
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
					if(c.sid != -1) {
						//cerr << c.sid;
					} else {
						//cerr << 'S';
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

	Dir get_dir(const Pos& from, const Pos& to) {
		//TODO en teoria no hace falta esta funcion si siempre son correctas las direcciones
		Dir d = to-from;
		for(int i = 0; i < all_dirs.size(); ++i) {
			if(d == all_dirs[i]) {
				return d;
			}
		}
		if(false) {
			return DEFAULT;
		}
		return {-1, -1};
	}

	bool dispara(const Starship &s) {
		Dir d;
		for(int i = 1; i <= 2; ++i) {
			d = {0, i};
			if(cell(s.pos+d).sid != -1 && player_of(cell(s.pos+d).sid) == me()) {
				if(i == 1) {
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

	void recalculate_route(const Starship &s, const Simulation &simulation) {
		Dir d = simulation.get_free_dir(s.pos, all_dirs);
		targets[s.sid].route = stack<Pos>();
		//cerr << "Se le da la direccion: " << d << endl;
		targets[s.sid].route.push(s.pos+d);
	}

	#define SIMULATE_STARSHIPS false

	void check_safe(const Starship &s) {
		int border = 2;
		Dir aux = {0, 2};
		if(!within_window(s.pos+aux, round())) {
			--border;
			aux = {0, 1};
			if(!within_window(s.pos+aux, round())) {
				--border;
			}
		}
		Simulation simulation = Simulation(s.pos, targets[s.sid].route.top(), number_rows(), m_universe, border);
		load_simulation(s.sid, simulation);
		//cerr << "inicio de la simulacion, estado inicio de ronda:" << endl;
		//print_simulation(simulation);
		simulation.run(SIMULATE_STARSHIPS, all_dirs);
		Cell c = simulation.get_cell(targets[s.sid].route.top());
		if(c.type == ASTEROID) {
			Cell aux_c;
			aux_c.type = EMPTY;
			simulation.set_cell(targets[s.sid].route.top(), aux_c);
		}
		//cerr << "fin de simulacion, resultado justo antes de mover la nave" << endl;
		//print_simulation(simulation);
		Dir d = get_dir(s.pos, targets[s.sid].route.top());
		aux = {-1, -1};
		if(d == aux || !simulation.can_move(s.pos, d)) {
			if(c.type == ASTEROID) {
				simulation.set_cell(s.pos+d, c);
			}
			//cerr << "Movimiento no seguro, recalcular ruta" << endl;
			recalculate_route(s, simulation);
		} else {
			//cerr << "Movimiento previsto seguro, no hay cambios" << endl;
		}
	}

	#define MIN_MISSILES 4

	void refresh_target(const Starship &s) {
		if(!targets[s.sid].ok || !is_target(targets[s.sid], targets[s.sid].n.first)) {
			TARGET_TYPE type = POINTS;
			if(s.nb_miss < MIN_MISSILES) {
				type = MISSILES;
				//cerr << "tiene que buscar misiles" << endl;
			}
			//TODO situarse mitad pantalla cuando esta muy adelante.
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
			//cerr << "No se le ha encontrdo objetivo, se le añade dir default" << endl;
		}
	}

	///
	/// \brief play    asigna un objetivo a cada nave al principio y mueve a cada nave seg?n su objetivo. También actualiza
	/// la posición de las naves enemigas en función de sus antiguas posiciones para no tener que hacer una b?squeda en cada turno.
	///

	virtual void play() {
		//cerr << "---------------------------" << endl;
		//cerr << "RONDA " << round() << endl;
		if(round() == 0) {
			m_universe = number_universe_columns();
			all_dirs = {FAST, FAST_UP, FAST_DOWN, DEFAULT, UP, DOWN, SLOW, SLOW_UP, SLOW_DOWN};
			targets = Targets(number_starships_per_player(), begin(me()), m_universe);
		}
		for(Starship_Id id = begin(me()); id != end(me()); ++id) {
			Starship s = starship(id);
			if(s.alive) {
				//cerr << "Turno de starship " << id << endl;
				//cerr << "Pos: ";
				//print_pos(s.pos);
				if(targets[id].route.empty()) {
					//cerr << "No tiene ruta" << endl;
					targets[id].ok = false;
				}

				refresh_target(s);
				check_safe(s);

				//cerr << "Su ruta actual es: " << endl;
				//print_route(s);

				if(cell(targets[id].route.top()).type == ASTEROID || cell(targets[id].route.top()).type == STARSHIP) {
					dispara(s);
				} else {
					move(id, get_dir(s.pos, targets[id].route.top()));
				}
				targets[id].next_pos = targets[id].route.top();
				targets[id].route.pop();
			} else {
				targets[id].ok = false;
			}
		}
		//TODO implementar update enemigos
	}
};

RegisterPlayer(PLAYER_NAME);

