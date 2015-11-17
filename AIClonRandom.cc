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

		Target operator[](Starship_Id id) {
			return v[id-offset];
		}

		void set_target(Starship_Id id, const Target &t) {
			v[id-offset] = t;
		}

	private:
		vector<Target> v;
		int offset = 0;
	};

	///
	/// \brief choose_target sirve para escoger un objetivo y un recorrido para
	/// \param s
	/// \return devuelve el objetivo seleccionado
	///
 /*
	int minimum_dist(matrix &map, pair<int, int> pos, int n, int m, int tesoros) {
		if(tesoros == 0) {
			return -1;
		}
		int dist = 0;
		queue<pair<int, int> > q;
		q.push(pos);
		map[pos.first][pos.second].second = true;
		int current_level = 1;
		int next_level = 0;
		while(!q.empty()) {
			pos = q.front();
			if(map[pos.first][pos.second].first == 't') {
				return dist;
			}
			q.pop();
			for(int i = pos.first-1; i <= pos.first+1; ++i) {
				if(i >= 0 && i < n) {
					for(int j = pos.second-(i == pos.first); j <= pos.second+(i == pos.first); ++j) {
						if(j >= 0 && j < m) {
							if(!map[i][j].second && map[i][j].first != 'X') {
								map[i][j].second = true;
								pair<int, int> pos_seg;
								pos_seg.first = i;
								pos_seg.second = j;
								q.push(pos_seg);
								++next_level;
							}
						}
					}
				}
			}
			--current_level;
			if(current_level == 0) {
				current_level = next_level;
				next_level = 0;
				++dist;
			}
		}
		return -1;
	}
	*/

	Target choose_target(const Starship &s) {
		Target t;
		int r = 0;
		queue<Pos> positions;
		positions.push(s.pos);
		set<Pos> visited;
		visited.insert(s.pos);
		while(positions.empty()) {

		}

		return t;
	}

	///
	/// \brief refresh_target comprueba que el objetivo sigue siendo accesible, existe y merece la pena, si no, escoge otro.
	/// Además, si en medio de un objetivo se da cuenta que puede conseguir algo de camino actualiza su recorrido.
	/// \param s
	///
	void refresh_target(const Starship &s) {
		if(cell(targets[s.sid].p).type == targets[s.sid].type) {
			targets.set_target(s.sid, choose_target(s));
		}
	}

	Targets targets;

	//TODO estructura de datos para almacenar naves enemigas

	///
	/// \brief play	asigna un objetivo a cada nave al principio y mueve a cada nave según su objetivo. También actualiza
	/// la posición de las naves enemigas en función de sus antiguas posiciones para no tener que hacer una búsqueda en cada turno.
	///
	virtual void play () {
		if(round() == 0) {
			targets = Targets(number_starships_per_player(), begin(me()));
			for(Starship_Id id = begin(me()); id != end(me()); ++id) {
				targets.set_target(id, choose_target(starship(id)));
			}
		}
		for(Starship_Id id = begin(me()); id != end(me()); ++id) {
			Starship s = starship(id);
			refresh_target(s);
			//TODO implementar refresh_target
		}
		//TODO implementar update enemigos
	}
};

RegisterPlayer(PLAYER_NAME);

