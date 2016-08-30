#include "RaptorApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"
#include "combust.h"
#include "heatconduction.h"
#include "heatconductiontimederivative.h"
#include "heatconvection.h"
#include "Reacted.h"
#include "Products.h"

template<>
InputParameters validParams<RaptorApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  params.set<bool>("use_legacy_output_syntax") = false;

  return params;
}

RaptorApp::RaptorApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  RaptorApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  RaptorApp::associateSyntax(_syntax, _action_factory);
}

RaptorApp::~RaptorApp()
{
}

// External entry point for dynamic application loading
extern "C" void RaptorApp__registerApps() { RaptorApp::registerApps(); }
void
RaptorApp::registerApps()
{
  registerApp(RaptorApp);
}

// External entry point for dynamic object registration
extern "C" void RaptorApp__registerObjects(Factory & factory) { RaptorApp::registerObjects(factory); }
void
RaptorApp::registerObjects(Factory & factory)
{

registerKernel(combust);
registerKernel(heatconduction);
registerKernel(heatconductiontimederivative);
registerKernel(heatconvection);
registerAuxKernel(Reacted);
registerMaterial(Products);
}

// External entry point for dynamic syntax association
extern "C" void RaptorApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { RaptorApp::associateSyntax(syntax, action_factory); }
void
RaptorApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
