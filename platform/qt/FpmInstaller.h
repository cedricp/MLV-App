#ifndef FPMINSTALLER_H
#define FPMINSTALLER_H

#include <QString>
#include <QStringList>
#include <QDir>
#include <QCoreApplication>
#include <QDebug>
#include <QMessageBox>

class FpmInstaller
{
public:
  static bool installFpm( QString fileName )
  {
      //Where to install it?
      QString newFileName = QCoreApplication::applicationDirPath().append( "/" ).append( QFileInfo(fileName).fileName() );
      //Remove existing script
      if( QFileInfo( newFileName ).exists() )
      {
          QFile( newFileName ).remove();
      }
      //Copy new one to app
      bool ret = QFile::copy( fileName, newFileName );
      return ret;
  }
  static void installFpm( QStringList *fileNameList )
  {
      for( int i = 0; i < fileNameList->count(); i++ )
      {
          QString fileName = fileNameList->at(i);
          //Where to install it?
          QString newFileName = QCoreApplication::applicationDirPath().append( "/" ).append( QFileInfo(fileName).fileName() );
          //Remove existing script
          if( QFileInfo( newFileName ).exists() )
          {
              QFile( newFileName ).remove();
          }
          //Copy new one to app
          bool ret = QFile::copy( fileName, newFileName );
          if( !ret ) fileNameList->removeAt(i);
      }
      return;
  }
};

#endif // FPMINSTALLER_H
