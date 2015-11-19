#include "Player.hh"
#include <vector>
#include <set>

using namespace std;

#define PLAYER_NAME ClonRandom

struct PLAYER_NAME : public Player {

	static Player* factory () {
		return new PLAYER_NAME;
	}

	///
	/// \brief The Target struct representa un objetivo escogido para una nave
	/// con su ruta optima, posicion final, su tipo y misiles necesarios.
	///
	struct Target {
		queue<Pos> ruta;
		Pos p;
		int misiles_nec;
		int rounds_nec;
		CType type;
	};

	///
	/// \brief The Targets class es un contenedor que almacena los objetivos de
	/// cada nave con alguna operación extra.
	///
	class Targets {
	public:
		Targets() {}

		Targets(int size, int offset) {
			this->offset = offset;
			v = vector<Target>(size);
		}

		Target& operator[](Starship_Id id) {
			return v[id-offset];
		}

	private:
		vector<Target> v;
		int offset = 0;
	};

	vector<Dir> all_dirs;

	struct compare {
		bool operator() (const Pos& p1, const Pos& p2) const{
			if(first(p1) != first(p2)) {
				return first(p1) < first(p2);
			} else {
				return second(p1) < second(p2);
			}
		}
	};

	bool cell_ok(const Pos &p) {
		Cell c = cell(p);
		return c.type == EMPTY || c.type == POINT_BONUS || c.type == MISSILE_BONUS;
	}

	int check_movement(const Pos &p, const Dir &d) {
		if(d == DEFAULT) {
			return true;
		} else if(!cell_ok(p+d)) {
			return false;
		} else if(d == FAST) {
			return cell_ok(p+DEFAULT);
		} else if(d == FAST_UP) {
			return check_movement(p, UP) && check_movement(p, FAST);
		} else if(d == FAST_DOWN) {
			return check_movement(p, DOWN) && check_movement(p, FAST);
		} else if(d == UP) {
			return cell_ok(p+UP) && cell_ok(p+DEFAULT);
		} else if(d == DOWN) {
			return cell_ok(p+DOWN) && cell_ok(p+DEFAULT);
		} else {
			return true;
		}
	}

	///
	/// \brief choose_target sirve para escoger un objetivo y un recorrido para
	/// \param s
	/// \return devuelve el objetivo seleccionado
	///

	#define MAX_LEVEL 15
	#define N_TARGETS 3
	bool scan_target(const Starship &s, Target &t, CType type) {
		vector<Target> v = vector<Target>(N_TARGETS);
		set<Pos, compare> visited = set<Pos, compare>();
		queue<Pos> q;
		q.push(s.pos);
		visited.insert(s.pos);
		int r = 0;
		int i = 0;
		int current_level = 1;
		int next_level = 0;
		Pos act;
		while(!q.empty() && i < N_TARGETS && r <= MAX_LEVEL) {
			act = q.front();
			if(cell(act).type == type) {
				v[i].p = act;
				v[i].rounds_nec = r;
				v[i].type = type;
				//v[i].type = cell(act).type;
				++i;
			}
			for(int j = 0; j < all_dirs.size()-1; ++j) {
				Pos new_pos = act+all_dirs[j];
				if(visited.find(new_pos) == visited.end()
						&& within_window(new_pos, round()+r)
						&& check_movement(act, all_dirs[j])) {
					visited.insert(new_pos);
					q.push(new_pos);
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
		if(i != 0) {
			t = v[0];
			return true;
		} else {
			return false;
		}
	}

	Target get_target(const Starship &s) {
		Target t;
		scan_target(s, t, POINT_BONUS);
		cerr << '(' << first(t.p) << " - " << second(t.p) << ')' << endl;
		cerr << t.type << endl;
		cerr << t.misiles_nec << endl;
		return t;
	}

	///
	/// \brief refresh_target comprueba que el objetivo sigue siendo accesible, existe y merece la pena, si no, escoge otro.
	/// Adem?s, si en medio de un objetivo se da cuenta que puede conseguir algo de camino actualiza su recorrido.
	/// \param s
	///
	void refresh_target(const Starship &s) {
		if(cell(targets[s.sid].p).type != targets[s.sid].type || !within_window(targets[s.sid].p, round())) {
			targets[s.sid] = get_target(s);
		}
	}

	Targets targets;

	//TODO estructura de datos para almacenar naves enemigas

	///
	/// \brief play    asigna un objetivo a cada nave al principio y mueve a cada nave seg?n su objetivo. También actualiza
	/// la posición de las naves enemigas en función de sus antiguas posiciones para no tener que hacer una b?squeda en cada turno.
	///
	virtual void play () {
		if(round() == 0) {
			all_dirs = {FAST, FAST_UP, FAST_DOWN, DEFAULT, UP, DOWN, SLOW_UP, SLOW_DOWN, SLOW};
			targets = Targets(number_starships_per_player(), begin(me()));
			for(Starship_Id id = begin(me()); id != end(me()); ++id) {
				targets[id] = get_target(starship(id));
				while(true);
			}
		}
		for(Starship_Id id = begin(me()); id != end(me()); ++id) {
			Starship s = starship(id);
			if(s.alive) {
				refresh_target(s);
				//TODO implementar refresh_target
			}
		}
		//TODO implementar update enemigos
	}
};

RegisterPlayer(PLAYER_NAME);
