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
	#define MAX_LEVEL_BFS 15
	#define N_TARGETS_BFS 4
	#define COST_MISSILE 0

	static Player* factory () {
		return new PLAYER_NAME;
	}

	int m_universe;

	const Dir INVALID_DIR = {-1, -1};

	struct Nodo_data {
		Pos prev;
		int nb_miss;
		int r;
		bool reduced;
	};

	typedef pair<Pos, Nodo_data> Nodo;

	struct compare_pos {
		bool operator() (const Pos& p1, const Pos& p2) const{
			if(first(p1) != first(p2)) {
				return first(p1) < first(p2);
			} else {
				return second(p1) < second(p2);
			}
		}
	};

	typedef map<Pos, Nodo_data, compare_pos> Nodos;

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
				//TODO si se llega antes con esta nave, se le quita el objetivo a la otra y se le asigna a esta.
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
										if(all_dirs[d] == DEFAULT || all_dirs[d] == FAST) {
											new_pos.second.sid = -2;
										} else {
											new_pos.second.sid = -3;
										}
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
			for(int d = 0; d < (int)all_dirs.size(); ++d) {
				pair<int, int> ij = get_ij(from+all_dirs[d]);
				if(affects(ij) && can_move(from, all_dirs[d], has_missiles && stage[ij.first][ij.second].sid > -2)) {
					return all_dirs[d];
				}
			}
			if(strict) {
				cerr << "No hay dir libre, se baja el nivel de strict" << endl;
				strict = false;
				Dir d = get_free_dir(from, all_dirs, slow, has_missiles);
				strict = true;
				return d;
			} else {
				cerr << "No hay ninguna dir libre, default" << endl;
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
				return c.type != ASTEROID && c.type != MISSILE && (c.type != STARSHIP || (!strict && c.sid <= -2));
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
		print_simulation(simulation);
		simulation.run(all_dirs);
		print_simulation(simulation);
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

	Pos choose_target(queue<Pos> &candidates) {
		//TODO escoger el target que tenga otro target del mismo tipo en menos de n (pequeña) num de rondas desde ahi.
		Pos p;
		int min_j = -1;
		while(!candidates.empty()) {
			if(min_j == -1 || second(candidates.front()) < min_j) {
				p = candidates.front();
				min_j = second(p);
			}
			candidates.pop();
		}
		return p;
	}

	bool is_target(const Target &t, const Pos &p, int r) {
		if(t.type == REPOSITION && (r == -1 || column_of_window(second(p), round()+r) <= column_of_window(round()+r+number_window_columns()*LIMIT_POS, round()+r))) {
			return true;
		} else if(t.type == POINTS && cell(p).type == POINT_BONUS) {
			return true;
		} else if(t.type == MISSILES && cell(p).type == MISSILE_BONUS) {
			return true;
		} else if(t.type == ANY && (cell(p).type == MISSILE_BONUS || cell(p).type == POINT_BONUS)) {
			return true;
		} else if(t.type == ENEMY && cell(p).type == STARSHIP) {
			return true;
		} else {
			return false;
		}
	}

	void scan_target(const Starship &s, const Target &t, queue<Pos> &candidates, Nodos &nodos) {
		queue<Pos> q = queue<Pos>();
		Nodo n;
		n.first = s.pos;
		n.second.prev = n.first;
		n.second.nb_miss = 0;
		n.second.r = -1;
		n.second.reduced = false;
		nodos.insert(n);
		q.push(s.pos);
		int r = 0;
		int current_level = 1;
		int next_level = 0;
		while(!q.empty() && candidates.size() < N_TARGETS_BFS && r <= MAX_LEVEL_BFS) {
			Pos act = q.front();
			q.pop();
			if(is_target(t, act, r) && targets.avaliable(s.sid, act)) {
				cerr << "Objetivo encontrado en ";
				print_pos(act);
				candidates.push(act);
			}
			for(int d = 0; d < (int)all_dirs.size(); ++d) {
				n.first = act+all_dirs[d];
				if(nodos.find(n.first) == nodos.end() && within_window(n.first, round()+r) && check_movement(act, all_dirs[d], s.nb_miss > 0)) {
					if(!SIMULATION_ON_SCAN || r != 0 || simulation.can_move(act, all_dirs[d], s.nb_miss > 0)) {
						nodos.insert(n);
						++next_level;
						q.push(n.first);
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
	}

	inline int cost(const Nodo &n) {
		return n.second.nb_miss*COST_MISSILE+n.second.r;
	}

	typedef pair<int, Nodos::iterator> arc;

	struct compare_arc {
		bool operator() (const arc& a, const arc& b) const{
			return a.first < b.first;
		}
	};

	void dijkstra(const Starship &s, Nodos::iterator i, Nodos &nodos) {
		i->second.r = 0;
		i->second.prev = i->first;
		priority_queue<arc, vector<arc>, compare_arc > q;
		q.push(arc(0, i));
		while(!q.empty()) {
			Nodos::iterator n = q.top().second;
			q.pop();
			if(!n->second.reduced) {
				n->second.reduced = true;
				for(int d = 0; d < (int)all_dirs.size(); ++d) {
					Nodos::iterator aux = nodos.find(n->first+all_dirs[d]);
					if(aux != nodos.end() && within_window(n->first, round()+n->second.r) && check_movement(n->first, all_dirs[d], s.nb_miss-n->second.nb_miss > 0)) {
						if(aux->second.r == -1 || cost(*aux) > cost(*n)+1+!cell_ok(aux->first)*COST_MISSILE) {
							if(round() == 79) {
								cerr << "para llegar a ";
								print_pos(aux->first);
								cerr << "se ha encontrado un coste menor desde ";
								print_pos(n->first);
							}
							aux->second.prev = n->first;
							aux->second.r = n->second.r+1;
							aux->second.nb_miss = n->second.nb_miss + !cell_ok(aux->first);
							q.push(arc(cost(*aux), aux));
						}
					}
				}
			}
		}
	}

	bool get_target(const Starship &s, Target &t) {
		t.reset_route();
		queue<Pos> candidates = queue<Pos>();
		Nodos nodos = Nodos();
		scan_target(s, t, candidates, nodos);
		if(candidates.empty()) {
			return false;
		} else {
			dijkstra(s, nodos.find(s.pos), nodos);
			t.n = *nodos.find(choose_target(candidates));
			t.set_route(nodos, t.n);
			return true;
		}
	}

	void print_pos(const Pos &p) {
		cerr << '(' << first(p) << ',' << second(p)%m_universe << ')' << endl;
	}

	void print_route(const Starship &s) {
		stack<Pos> pila = targets[s.sid].route;
		while(!pila.empty()) {
			print_pos(pila.top());
			pila.pop();
		}
	}

	void print_simulation(const Simulation &s) {
		for(int i = 0; i < (int)s.stage.size(); ++i) {
			for(int j = 0; j < (int)s.stage[0].size(); ++j) {
				Cell c = s.stage[i][j];
				if(c.type == ASTEROID) {
					cerr << 'X';
				} else if(c.type == POINT_BONUS) {
					cerr << 'P';
				} else if(c.type == MISSILE_BONUS) {
					cerr << 'B';
				} else if(c.type == MISSILE) {
					cerr << 'M';
				} else if(c.type == STARSHIP) {
					if(c.sid == -1) {
						cerr << 'S';
					} else if(c.sid <= -2) {
						cerr << 's';
					} else {
						cerr << c.sid;
					}
				} else {
					cerr << '.';
				}
			}
			cerr << endl;
		}
		cerr << endl;
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
		targets[s.sid].ok = get_target(s, targets[s.sid]);
		if(!targets[s.sid].ok) {
			cerr << "No se ha encontrado ningun objetivo accesible" << endl;
			targets[s.sid].reset_route();
			Dir d = simulation.get_free_dir(s.pos, all_dirs, targets[s.sid].type == REPOSITION, s.nb_miss > 0);
			cerr << "Direccion para escapar: " << d << endl;
			targets[s.sid].route.push(s.pos+d);
			targets[s.sid].ok = true;
		}
	}

	void check_safe(const Starship &s) {
		Dir d = get_dir(s.pos, targets[s.sid].route.top());
		if(!simulation.can_move(s.pos, d, s.nb_miss > 0)) {
			cerr << "Movimiento no seguro, recalcular ruta" << endl;
			recalculate_route(s);
		} else {
			cerr << "Movimiento previsto seguro, no hay cambios" << endl;
		}
	}

	void refresh_target(const Starship &s) {
		if(!targets[s.sid].ok || !is_target(targets[s.sid], targets[s.sid].n.first, -1) || get_dir(s.pos, targets[s.sid].route.top()) == INVALID_DIR || s.nb_miss < targets[s.sid].n.second.nb_miss) {
			TARGET_TYPE type = POINTS;
			if(s.nb_miss < MIN_MISSILES) {
				type = MISSILES;
			}
			if(column_of_window(second(s.pos), round()) >= column_of_window(round()+number_window_columns()*LIMIT_POS, round())) {
				cerr << "demasiado adelante, objetivo: reposicion" << endl;
				//type = REPOSITION;
			}
			cerr << "Buscar objetivo: " << type << endl;
			targets[s.sid].type = type;
			targets[s.sid].ok = get_target(s, targets[s.sid]);
			if(!targets[s.sid].ok) {
				targets[s.sid].type = ANY;
				targets[s.sid].ok = get_target(s, targets[s.sid]);
			} else {
				cerr << "Objetivo encontrado en: ";
				print_pos(targets[s.sid].n.first);
			}
		}
		if(!targets[s.sid].ok) {
			targets[s.sid].route.push(s.pos+DEFAULT);
			targets[s.sid].ok = true;
			cerr << "No se le ha encontrdo objetivo, se le añade dir SLOW" << endl;
		}
	}

	///
	/// \brief play    asigna un objetivo a cada nave al principio y mueve a cada nave seg?n su objetivo. También actualiza
	/// la posición de las naves enemigas en función de sus antiguas posiciones para no tener que hacer una b?squeda en cada turno.
	///

	virtual void play() {
		cerr << "--------------------------------------------" << endl;
		cerr << "RONDA " << round() << endl;
		if(round() == 0) {
			m_universe = number_universe_columns();
			all_dirs = {FAST, FAST_UP, FAST_DOWN, DEFAULT, UP, DOWN, SLOW_UP, SLOW_DOWN, SLOW};
			targets = Targets(number_starships_per_player(), begin(me()), m_universe);
		}
		for(Starship_Id id = begin(me()); id != end(me()); ++id) {
			const Starship s = starship(id);
			if(s.alive) {
				cerr << "--------------------------------------------" << endl;
				cerr << "Turno de starship " << id << endl;
				cerr << "Pos: ";
				print_pos(s.pos);
				if(targets[id].route.empty()) {
					cerr << "No tiene ruta" << endl;
					targets[id].ok = false;
				}

				load_simulation(s);
				refresh_target(s);
				check_safe(s);

				cerr << "Su ruta actual es: " << endl;
				print_route(s);
				cerr << "Dir: " << get_dir(s.pos, targets[id].route.top()) << endl;
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

