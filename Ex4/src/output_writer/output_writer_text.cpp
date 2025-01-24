#include "output_writer/output_writer_text.h"

void OutputWriterText::writeFile(double currentTime)
{
  // Assemble the filename
  std::stringstream fileName;
  fileName << "out/output/output_" << std::setw(4) << std::setfill('0') << fileNo_ << '.' 
          << rank_ << ".txt";
  
  // increment file no.
  fileNo_++;

  // open file
  std::ofstream file(fileName.str().c_str());

  if (!file.is_open())
  {
    std::cout << "Could not write to file \"" << fileName.str() << "\".";
    return;
  }

  // write time
  file << "t: " << currentTime << std::endl;

  const double dx = discretization_->meshWidth()[0];
  const double dy = discretization_->meshWidth()[1];
  const int nx = discretization_->nCells()[0];
  const int ny = discretization_->nCells()[1];
  // write mesh width
  file << "nCells: " << nx << "x" << ny << ", dx: " << dx << ", dy: " << dy << std::endl << std::endl;

  const int fieldWidth = 9;   // number of characters to use for a single value

  // write u
  // ---------
  // write header lines
  file << "u (" << discretization_->u().size()[0] << "x" << discretization_->u().size()[1] << "): " << std::endl 
    << std::string(fieldWidth, ' ') << "|";
  for (int i = -discretization_->ui0(); i < discretization_->uiN()+1; i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl << std::string(fieldWidth*(discretization_->u().size()[0]+2)+1, '-') << std::endl;

  // write u values
  for (int j = discretization_->ujN(); j >= -discretization_->uj0(); j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = -discretization_->ui0(); i < discretization_->uiN()+1; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << discretization_->u(i,j);
    }
    file << std::endl;
  }
  file << std::endl;

  // write v
  // ---------
  // write header lines
  file << "v (" << discretization_->v().size()[0] << "x" << discretization_->v().size()[1] << "): " << std::endl 
    << std::string(fieldWidth, ' ') << "|";
  for (int i = -discretization_->vi0(); i < discretization_->viN()+1; i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl << std::string(fieldWidth*(discretization_->v().size()[0]+2)+1, '-') << std::endl;

  // write v values
  for (int j = discretization_->vjN(); j >= -discretization_->vj0(); j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = -discretization_->vi0(); i < discretization_->viN()+1; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << discretization_->v(i,j);
    }
    file << std::endl;
  }
  file << std::endl;

  // write p
  // ---------
  // write header lines
  file << "p (" << discretization_->p().size()[0] << "x" << discretization_->p().size()[1] << "): " << std::endl 
    << std::string(fieldWidth, ' ') << "|";
  for (int i = -discretization_->pi0(); i < discretization_->piN()+1; i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl << std::string(fieldWidth*(discretization_->p().size()[0]+2)+1, '-') << std::endl;

  // write p values
  for (int j = discretization_->pjN(); j >= -discretization_->pi0(); j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = -discretization_->pi0(); i < discretization_->piN()+1; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << discretization_->p(i,j);
    }
    file << std::endl;
  }
  file << std::endl;

  // write f
  // ---------
  // write header lines
  file << "F (" << discretization_->u().size()[0] << "x" << discretization_->u().size()[1] << "): " << std::endl 
    << std::string(fieldWidth, ' ') << "|";
  for (int i = -discretization_->ui0(); i < discretization_->uiN()+1; i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl << std::string(fieldWidth*(discretization_->u().size()[0]+2)+1, '-') << std::endl;

  // write f values
  for (int j = discretization_->ujN(); j >= -discretization_->uj0(); j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = -discretization_->ui0(); i < discretization_->uiN()+1; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << discretization_->f(i,j);
    }
    file << std::endl;
  }
  file << std::endl;

  // write g
  // ---------
  // write header lines
  file << "G (" << discretization_->v().size()[0] << "x" << discretization_->v().size()[1] << "): " << std::endl 
    << std::string(fieldWidth, ' ') << "|";
  for (int i = -discretization_->vi0(); i < discretization_->viN()+1; i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl << std::string(fieldWidth*(discretization_->v().size()[0]+2)+1, '-') << std::endl;

  // write g values
  for (int j = discretization_->vjN(); j >= -discretization_->vj0(); j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = -discretization_->vi0(); i < discretization_->viN()+1; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << discretization_->g(i,j);
    }
    file << std::endl;
  }
  file << std::endl;

  // write rhs
  // ---------
  // write header lines
  file << "rhs (" << discretization_->rhs().size()[0] << "x" << discretization_->rhs().size()[1] << "): " << std::endl 
    << std::string(fieldWidth, ' ') << "|";
  // ???
  for (int i = 0; i < discretization_->piN(); i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl << std::string(fieldWidth*(discretization_->p().size()[0]+2)+1, '-') << std::endl;

  // write rhs values
  for (int j = discretization_->pjN()-1; j >= 0; j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = 0; i < discretization_->piN(); i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << discretization_->rhs(i,j);
    }
    file << std::endl;
  }
  file << std::endl;

}

// wellp, the numMaxPressureSteps is a const, so the clean way would be to define it on object definition
// but we use two seperate debuggers for vel and p, so two seperate output-writer types would then be the way to go
// this is very much a ductape way to implement a nice output file, which should suffice, given it is only for debugging
void OutputWriterText::writePressureFile(int simStep, int pressureStep, int numMaxPressureSteps, char substep)
{
  // counter for files, counter value is part of the file name
  static int pressurefileNo = 0;

  const int pressureNoExp = std::ceil(std::log10(numMaxPressureSteps));

  // Assemble the filename
  std::stringstream fileName;
  //fileName << "out/pressure_" << std::setw(4) << std::setfill('0') << pressurefileNo++ << ".txt";
  fileName << "out/pressure/pressure_" << std::setw(4) << std::setfill('0') << simStep << '.' 
           << std::setw(pressureNoExp) << std::setfill('0') << pressureStep << substep << '.' << rank_ << ".txt";
  
  // open file
  std::ofstream file(fileName.str().c_str());
  int returnValue = system("mkdir -p out/output");
  if (!file.is_open())
  {
    std::cout << "Could not write to file \"" << fileName.str() << "\".\n";
    return;
  }

  const double dx = discretization_->meshWidth()[0];
  const double dy = discretization_->meshWidth()[1];
  const int nx = discretization_->nCells()[0];
  const int ny = discretization_->nCells()[1];

  // write mesh width
  file << "nCells: " << discretization_->nCells()[0] << "x" << discretization_->nCells()[1] 
    << ", dx: " << dx << ", dy: " << dy << std::endl << std::endl;

  const int fieldWidth = 9;   // number of characters to use for a single value

  // write p
  // ---------
  // write header lines
  file << "p (" << discretization_->p().size()[0] << "x" << discretization_->p().size()[1] << "): " << std::endl 
    << std::string(fieldWidth, ' ') << "|";
  for (int i = -discretization_->pi0(); i < discretization_->piN()+1; i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl << std::string(fieldWidth*(discretization_->p().size()[0]+2)+1, '-') << std::endl;

  // write p values
  for (int j = discretization_->pjN(); j >= -discretization_->pj0(); j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = -discretization_->pi0(); i < discretization_->piN()+1; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << discretization_->p(i,j);
    }
    file << std::endl;
  }
  file << std::endl;

}
