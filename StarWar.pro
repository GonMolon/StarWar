TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    Action.cc \
    AIDemo.cc \
    AINull.cc \
    Board.cc \
    Game.cc \
    Main.cc \
    Player.cc \
    Registry.cc \
    Utils.cc \
    AIEsUnaTrampa.cc

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    Action.hh \
    Board.hh \
    Defs.hh \
    Game.hh \
    Player.hh \
    Registry.hh \
    Utils.hh

