#include "Player.hh"
#include <vector>
#include <map>
#include <stack>

using namespace std;

#define PLAYER_NAME ClonRandom

struct PLAYER_NAME : public Player {

	#define MIN_MISSILES 4
	#define STRICT_SIMULATION true
	#define BIG_SIMULATION true
	#define SIMULATION_ON_SCAN true
	#define LIMIT_POS 0.9 //[0-1] 1 = no limit
	#define MAX_LEVEL_BFS 30
	#define N_TARGETS_BFS 1

	static Player* factory () {
		return new PLAYER_NAME;
	}

	int m_universe;

	const Dir INVALID_DIR = {-1, -1};

	struct Data_pos {
		Pos prev;
		int nb_miss;
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
		bool ok;
		TARGET_TYPE type;
		Pos next_pos;

		void set_route(Nodos &nodos, Nodo n) {
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

		bool avaliable(Starship_Id id, Pos p) {
			p = {first(p), second(p)%m_universe};
			for(int i = 0; i < (int)v.size(); ++i) {
				if(id-offset != i && v[i].ok && v[i].type != ENEMY && first(v[i].n.first) == first(p) && second(v[i].n.first)%m_universe == second(p)%m_universe) {
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
				stack<pair<pair<int, int>, Cell> > cells = stack<pair<pair<int, int>, Cell> >();
				for(int j = (int)stage[0].size()-1; j >= 0; --j) {
					for(int i = 0; i < (int)stage.size(); ++i) {
						if(stage[i][j].type == STARSHIP && stage[i][j].sid == -1) {
							for(int d = 0; d < (int)all_dirs.size(); ++d) {
								pair<int, int> ij = {i, j};
								if(check_movement(ij, all_dirs[d], stage[i][j].mid > 0)) {
									ij.first += first(all_dirs[d]);
									ij.second += second(all_dirs[d]);
									if(affects(ij)) {
										pair<pair<int, int>, Cell> new_pos;
										new_pos.first = ij;
										new_pos.second = stage[i][j];
										if(new_pos.second.type == ASTEROID) {
											--new_pos.second.mid;
										}
										if(all_dirs[d] == DEFAULT || all_dirs[d] == FAST) {
											new_pos.second.sid = -2;
										} else {
											new_pos.second.sid = -3;
										}
										cells.push(new_pos);
									}
								}
							}

						}
					}
				}
				while(!cells.empty()) {
					stage[cells.top().first.first][cells.top().first.second] = cells.top().second;
					cells.pop();
				}
			}
		}

		bool can_move(Pos p, Dir d, bool has_missiles) const {
			reduce_pos(p);
			pair<int, int> ij = get_ij(p);
			if(!check_movement(ij, d, has_missiles && stage[get_ij(p+d).first][get_ij(p+d).second].sid > -2)) {
				return false;
			}
			p = p+d;
			for(int j = -2; j < 0; ++j) {
				Dir aux = {0, j};
				pair<int, int> ij = get_ij(p+aux);
				if(affects(ij)) {
					Cell c = stage[ij.first][ij.second];
					if(c.type == MISSILE || (c.type == STARSHIP && (c.sid == -1 || c.sid == -2 || (c.sid == -3 && strict)) && c.mid > 0)) {
						if(j == -1 || (stage[ij.first][ij.second+1].type == EMPTY || (d != SLOW_UP && d != SLOW_DOWN)) || stage[ij.first][ij.second+1].type == STARSHIP) {
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
				//cerr << "probando con " << all_dirs[d] << endl;
				pair<int, int> ij = get_ij(from+all_dirs[d]);
				if(affects(ij) && can_move(from, all_dirs[d], has_missiles && stage[ij.first][ij.second].sid > -2)) {
					//cerr << "VALE" << endl;
					return all_dirs[d];
				}
				//cerr << "no vale" << endl;
			}
			if(strict) {
				//cerr << "No hay dir libre, se baja el nivel de strict" << endl;
				strict = false;
				Dir d = get_free_dir(from, all_dirs, slow, has_missiles);
				strict = true;
				return d;
			} else {
				//cerr << "No hay ninguna dir libre, default" << endl;
				return DEFAULT;
			}
		}

		void add_starship(Starship_Id id, Pos p) {
			reduce_pos(p);
			pair<int, int> ij = get_ij(p);
			if(affects(ij) && stage[ij.first][ij.second].type != MISSILE) {
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
			Starship_Id id = stage[ij.first][ij.second].sid;
			ij.first += first(d);
			ij.second += second(d);
			if(affects(ij)) {
				Cell c = stage[ij.first][ij.second];
				return c.type != ASTEROID && c.type != MISSILE && (c.type != STARSHIP || (id >= 0 && (!strict && c.sid <= -2)) || (id < 0 && c.sid <= -2));
			} else {
				return false;
			}
		}

		bool check_movement(const pair<int, int> &ij, const Dir &d, bool has_missiles) const {
			if(d == DEFAULT) {
				return cell_ok(ij, d) || (has_missiles && stage[ij.first][ij.second+1].sid < 0);
			} else if(d != SLOW && !cell_ok(ij, d)) {
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
		if(!within_window(s.pos+aux, round()+1)) {
			--border;
			aux = {0, 1};
			if(!within_window(s.pos+aux, round()+1)) {
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
		for(Starship_Id id_aux = begin(me()); id_aux != s.sid; ++id_aux) {
			if(id_aux != s.sid && targets[id_aux].ok) {
				simulation.add_starship(id_aux, targets[id_aux].next_pos);
			}
		}
		//print_simulation(simulation);
		simulation.run(all_dirs);
		//print_simulation(simulation);
	}

	vector<Dir> all_dirs;

	inline int column_of_window(int j, int r) {
			j = j%m_universe;
			int aux = r%m_universe;
			if(j < aux) {
					j += m_universe;
			}
			return j-aux;
	}

	bool cell_ok(const Pos &p) {
		return cell(p).type != ASTEROID;
	}

	bool check_movement(const Pos &p, const Dir &d, bool has_missiles) {
		if(d == DEFAULT) {
			return cell_ok(p+d) || has_missiles;
		} else if(d != SLOW && !cell_ok(p+d)) {
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

	Nodo choose_target(queue<Nodo> &candidates) {
		Nodo n = candidates.front();
		candidates.pop();
		while(!candidates.empty()) {
			if(second(candidates.front().first) < second(n.first)) {
				n = candidates.front();
			}
			candidates.pop();
		}
		return n;
	}

	bool is_target(const Target &t, const Pos &p, int r) {
		CType type = cell(p).type;
		if(t.type == REPOSITION && (r == -1 || (column_of_window(second(p), r) <= column_of_window(r+number_window_columns()*LIMIT_POS, r) || type == POINT_BONUS || type == MISSILE_BONUS))) {
			return true;
		} else if(t.type == POINTS && type == POINT_BONUS) {
			return true;
		} else if(t.type == MISSILES && type == MISSILE_BONUS) {
			return true;
		} else if(t.type == ANY && (type == MISSILE_BONUS || type == POINT_BONUS)) {
			return true;
		} else if(t.type == ENEMY && type == STARSHIP and is_enemy(p)) {
			return true;
		}
		return false;
	}

	static int priorize_dir(const vector<Dir> &all_dirs, int d, bool slow) {
		if(!slow) {
			return d;
		} else {
			if(all_dirs[d] == FAST || all_dirs[d] == FAST_UP || all_dirs[d] == FAST_DOWN) {
				return all_dirs.size()-1;
			}
			return all_dirs.size()-1-d;
		}
	}

	//TODO realizar adaptacion alternativa algoritmo dijkstra: busqueda+camino minimo a la vez
	Nodo scan_target(const Starship &s, Target &t, Nodos &visited) {
		bool slow = t.type == REPOSITION;
		queue<Nodo> q;
		queue<Nodo> candidates;
		Nodo next;
		next.first = s.pos;
		next.second.prev = next.first;
		next.second.nb_miss = 0;
		q.push(next);
		visited.insert(next);
		int r = 0;
		int current_level = 1;
		int next_level = 0;
		while(!q.empty() && candidates.size() < N_TARGETS_BFS && r <= MAX_LEVEL_BFS) {
			Nodo act = q.front();
			q.pop();
			if(is_target(t, act.first, r) && targets.avaliable(s.sid, act.first)) {
				candidates.push(act);
			} else {
				for(int j = 0; j < (int)all_dirs.size(); ++j) {
					int d = priorize_dir(all_dirs, j, slow);
					next.first = act.first+all_dirs[d];
					if((visited.find(next.first) == visited.end())
							&& within_window(next.first, round()+r)
							&& check_movement(act.first, all_dirs[d], s.nb_miss-act.second.nb_miss > 0)) {
						if(!SIMULATION_ON_SCAN || r != 0 || simulation.can_move(act.first, all_dirs[d], s.nb_miss > 0)) {
							next.second.prev = act.first;
							next.second.nb_miss = act.second.nb_miss + !cell_ok(next.first);
							visited.insert(next);
							q.push(next);
							++next_level;
						}
					}
				}
			}
			--current_level;
			if(current_level == 0) {
				current_level = next_level;
				next_level = 0;
				++r;
			}
		}
		if(!candidates.empty()) {
			return choose_target(candidates);
		} else {
			Nodo aux;
			aux.first = {-1, -1};
			return aux;
		}
	}

	bool get_target(const Starship &s, Target &t) {
		t.reset_route();
		Nodos visited = Nodos();
		Nodo n = scan_target(s, t, visited);
		Pos not_found = {-1, -1};
		if(n.first == not_found) {
			t.ok = false;
			return false;
		} else {
			t.set_route(visited, n);
			if(t.route.size() == 0) {
				t.route.push(s.pos);
			}
			n.first = {first(n.first), second(n.first)%m_universe};
			t.n = n;
			t.ok = true;
			return true;
		}
	}

	void print_pos(const Pos &p) {
		//cerr << '(' << first(p) << ',' << second(p)%m_universe << ')' << endl;
	}

	void print_route(const Starship &s) {
		stack<Pos> pila = targets[s.sid].route;
		while(!pila.empty()) {
			////print_pos(pila.top());
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
					if(c.sid == -1) {
						//cerr << 'S';
					} else if(c.sid <= -2) {
						//cerr << 's';
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

	inline bool is_enemy(const Pos &p) {
		Cell c = cell(p);
		if(c.type != STARSHIP) {
			return false;
		} else {
			return player_of(c.sid) != me();
		}
	}

	inline bool shoot_enemy(const Starship &s) {
		return s.nb_miss > 0 && (is_enemy(s.pos+DEFAULT) || (cell(s.pos+DEFAULT).type == EMPTY && is_enemy(s.pos+FAST)));
	}

	void recalculate_route(const Starship &s) {
		targets[s.sid].type = ANY;
		get_target(s, targets[s.sid]);
		if(!targets[s.sid].ok) {
			//cerr << "No se ha encontrado ningun objetivo accesible" << endl;
			targets[s.sid].reset_route();
			Dir d = simulation.get_free_dir(s.pos, all_dirs, targets[s.sid].type == REPOSITION, s.nb_miss > 0);
			//cerr << "Direccion para escapar: " << d << endl;
			targets[s.sid].route.push(s.pos+d);
			targets[s.sid].n.first = targets[s.sid].route.top();
			targets[s.sid].ok = true;
		}
	}

	void check_safe(const Starship &s) {
		Dir d = get_dir(s.pos, targets[s.sid].route.top());
		if(!simulation.can_move(s.pos, d, s.nb_miss > 0)) {
			//cerr << "Movimiento no seguro, recalcular ruta" << endl;
			recalculate_route(s);
		}
	}

	void refresh_target(const Starship &s) {
		if(!targets[s.sid].ok || !is_target(targets[s.sid], targets[s.sid].n.first, -1) || get_dir(s.pos, targets[s.sid].route.top()) == INVALID_DIR || s.nb_miss < targets[s.sid].missiles_nec) {
			TARGET_TYPE type;
			type = POINTS;
			if(s.nb_miss < MIN_MISSILES) {
				type = MISSILES;
			}

			if(column_of_window(second(s.pos), round()) > column_of_window(round()+number_window_columns()*LIMIT_POS, round())) {
					//type = REPOSITION;
			}
			targets[s.sid].type = type;
			get_target(s, targets[s.sid]);
			if(!targets[s.sid].ok) {
				targets[s.sid].type = ANY;
				get_target(s, targets[s.sid]);
				if(!targets[s.sid].ok) {
					targets[s.sid].reset_route();
					if(type == REPOSITION) {
						targets[s.sid].route.push(s.pos+SLOW);
					} else {
						targets[s.sid].route.push(s.pos+DEFAULT);
					}
					targets[s.sid].n.first = targets[s.sid].route.top();
					targets[s.sid].ok = true;
				}
			}
		}
	}

	bool save_cpu = false;

	virtual void play() {
		//cerr << "-----------------------------" << endl;
		//cerr << "RONDA " << round() << endl;
		//cerr << "-----------------------------" << endl;
		if(round() == 0) {
			m_universe = number_universe_columns();
			all_dirs = {FAST, FAST_UP, FAST_DOWN, DEFAULT, UP, DOWN, SLOW_UP, SLOW_DOWN, SLOW};
			targets = Targets(number_starships_per_player(), begin(me()), m_universe);
		} else {
			if(!save_cpu && status(me())*number_rounds()/round() > 0.9) {
				save_cpu = true;
			}
		}

		for(Starship_Id id = begin(me()); id != end(me()); ++id) {
			const Starship s = starship(id);
			if(s.alive) {
				//cerr << "Turno de starship " << id << endl;
				//cerr << "Pos: ";
				//print_pos(s.pos);

				if(!save_cpu) {
					targets[id].ok = false; //mejores resultados pero mas cpu. solamente usar si sobra
				}
				if(targets[id].route.empty()) {
					//cerr << "No tiene ruta" << endl;
					targets[id].ok = false;
				}

				load_simulation(s);
				refresh_target(s);
				check_safe(s);
				//cerr << "Su ruta actual es: " << endl;
				//print_route(s);
				//cerr << "Dir: " << get_dir(s.pos, targets[id].route.top()) << endl;
				if(shoot_enemy(s) || cell(targets[id].route.top()).type == ASTEROID) {
					shoot(id);
					targets[s.sid].next_pos = s.pos+DEFAULT;
				} else {
					move(id, get_dir(s.pos, targets[id].route.top()));
					targets[id].next_pos = targets[id].route.top();
				}
				targets[id].route.pop();
			} else {
				targets[id].ok = false;
			}
		}
	}
};

///
/// TODO Tengo que poder permanecer en la misma casilla al hacer una ruta: MEJORAR BFS. DIJKSTRA
/// TODO Mejorar posicionamiento en ventana
/// TODO Asignar estrategia a cada nave: Una recolectora y otra cazadora que proteja a la primera. Activar cazadora despues de 30% de partida
/// TODO Objetivo starship: Predecir a que objetivo va a ir mediante un bfs simple e ir hacia alli
/// TODO Decidir objetivo en funcion de los que se puedan coger desde ahi (bfs dentro de bfs)
/// TODO Calcular valor medio + desviacion tipica de la j de las posiciones enemigas para determinar posicion de pantalla.
/// TODO Mejorar esquive de misiles en situaciones extremas con muchas naves rodeando. (Partida Sanfe).
///

RegisterPlayer(PLAYER_NAME);
