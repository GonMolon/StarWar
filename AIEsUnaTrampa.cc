#include "Player.hh"
#include <vector>

using namespace std;

#define PLAYER_NAME EsUnaTrampa

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
	Target choose_target(const Starship &s) {
		Target t;
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
		}
	}
};

RegisterPlayer(PLAYER_NAME);

