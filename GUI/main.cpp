/****************************************************************************
****************************************************************************/

//! [0]
#include <QApplication>

#include "mainwindow.h"
#include "log.h"
#include <windows.h>
#include <stdio.h>

LOG_DECLARE;

/*
BOOL CtrlHandler( DWORD fdwCtrlType )
{
  char msg[1024];

  switch( fdwCtrlType )
  {
    // Handle the CTRL-C signal.
    case CTRL_C_EVENT:
      sprintf(msg, "Ctrl-C event\n\n" );
      LOG_MSG(msg);
      Beep( 750, 300 );
      return( TRUE );

    // CTRL-CLOSE: confirm that the user wants to exit.
    case CTRL_CLOSE_EVENT:
      Beep( 600, 200 );
      sprintf(msg, "Ctrl-Close event\n\n" );
      LOG_MSG(msg);
      return( TRUE );

    // Pass other signals to the next handler.
    case CTRL_BREAK_EVENT:
      Beep( 900, 200 );
      sprintf(msg, "Ctrl-Break event\n\n" );
      LOG_MSG(msg);
      return FALSE;

    case CTRL_LOGOFF_EVENT:
      Beep( 1000, 200 );
      sprintf(msg, "Ctrl-Logoff event\n\n" );
      LOG_MSG(msg);
      return FALSE;

    case CTRL_SHUTDOWN_EVENT:
      Beep( 750, 500 );
      sprintf(msg, "Ctrl-Shutdown event\n\n" );
      LOG_MSG(msg);
      return FALSE;

    default:
      LOG_MSG("CtrlHandler got an event:"); // %d",fdwCtrlType)
      return FALSE;
  }
}
*/

int main(int argc, char *argv[])
{
    int res;
    char msg[1024];
    //initialize file logger
    LOG_INIT("spheroid_GUI.log");

    /*
    if( SetConsoleCtrlHandler( (PHANDLER_ROUTINE) CtrlHandler, TRUE ) )
    {
        sprintf(msg, "\nThe Control Handler is installed.\n" );
        LOG_MSG(msg);
//      printf( "\nThe Control Handler is installed.\n" );
//      printf( "\n -- Now try pressing Ctrl+C or Ctrl+Break, or" );
//      printf( "\n    try logging off or closing the console...\n" );
//      printf( "\n(...waiting in a loop for events...)\n\n" );

//      while( 1 ){ }
    }
    else
    {
      sprintf(msg, "\nERROR: Could not set control handler");
      LOG_MSG(msg);
      return 1;
    }
    */

    QApplication app(argc, argv);

    MainWindow mainWin;
    mainWin.show();

    res = app.exec();
    sprintf(msg,"Result code: %d",res);
    LOG_MSG(msg);
    return res;
}
//! [0]
