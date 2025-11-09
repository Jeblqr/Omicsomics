import React from 'react';
import { ProjectSwitcher } from '../../components/ProjectSwitcher';
import { SandboxView } from '../../components/SandboxView';

export const SandboxPage: React.FC = () => {
  return (
    <div>
      <ProjectSwitcher />
      <SandboxView />
    </div>
  );
};
